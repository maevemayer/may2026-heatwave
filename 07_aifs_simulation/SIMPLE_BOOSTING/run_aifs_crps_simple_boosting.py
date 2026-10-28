import os
import pickle
import glob
import argparse
from anemoi.inference.runners.simple import SimpleRunner


def run_aifs(
    input_path,
    output_path,
    date_str,
    exp_name,
    checkpoint,
    lead_time=720,
    target_variables=["2t", "10u", "10v", "z_500"],
    times_to_save=[0, 6, 12, 18],
    N_MEMBERS=50,
    mask=None,
):
    """Runs the AIFS model for a given set of initial states and saves the forecasts.

    Args:
        input_path (str): Directory containing the input .pkl files with initial states.
        output_path (str): Directory where the output .pkl files with forecasts will be saved.
        date_str (str): Date string in the format "YYYYMMDDHH" corresponding to the initial states.
        checkpoint (str): Path to the AIFS model checkpoint file.
        lead_time (int, optional): Lead time in hours for the forecast. Defaults to 720 (30 days).
        target_variables (list[str], optional): List of variable names to include in the output forecasts. Defaults to ["2t", "10u", "10v", "z_500"].
        times_to_save (list[int], optional): List of hours (0-23) at which to save the forecast steps. Defaults to [0, 6, 12, 18] (every 6 hours).
        N_MEMBERS (int, optional): Number of ensemble members to process (not including the control member 00). Defaults to 50
        mask (np.array[boolean]): a mask indicating on which zone should the outputs be saved

    Returns:
        None, but saves the forecast outputs as .pkl files in the specified output directory.
    """
    # load files
    dir_name=f"{date_str}-{exp_name}"
    dir_pattern = os.path.join(input_path, f"{date_str}-{exp_name}")
    found_dir = glob.glob(dir_pattern)[0]
    if not found_dir:
        raise FileNotFoundError(f"No input directory found matching: {dir_pattern}")

    # creates the output directry if it doesn't exist
    out_dir = os.path.join(output_path, dir_name)
    os.makedirs(out_dir, exist_ok=True)

    input_filename = os.path.join(found_dir, f"init_state-{dir_name}.pkl")
    # load the initial state from the input file
    with open(input_filename, "rb") as f:
        input_state = pickle.load(f)

    # loading the model
    runner = SimpleRunner(checkpoint, device="cuda")

    # for each member, load the initial state, run the model and save the output
    for i in range(N_MEMBERS + 1):
        # format the member id to be two digits, e.g. 00, 01, ..., 50 and construct the input filename
        member_id = f"{i:02d}"
        print(f"Running simulation for member: {member_id}",flush=True)
        
        # run the model and collect the forecast steps
        forecast_steps = []
        # save initial state as the first step of the forecast
        forecast_steps.append(
            {
                "date": input_state["date"],
                "fields": {
                    k: v[1][mask] if mask is not None else v[1]
                    for k, v in input_state["fields"].items()
                    if k in target_variables
                },
            }
        )
        # forecast
        for state in runner.run(input_state=input_state, lead_time=lead_time):
            # filter the state to only include the target variables
            if state["date"].hour in times_to_save:
                filtered_state = {
                    "date": state["date"],
                    "fields": {
                        k: v[mask] if mask is not None else v
                        for k, v in state["fields"].items()
                        if k in target_variables
                    },
                }
                forecast_steps.append(filtered_state)

        # save the forecast steps to the output file
        output_filename = os.path.join(
            out_dir, f"simple_boosting-{dir_name}-{member_id}.pkl"
        )
        with open(output_filename, "wb") as f_out:
            pickle.dump(forecast_steps, f_out)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_path", type=str, required=True)
    parser.add_argument("--output_path", type=str, required=True)
    parser.add_argument("--date_str", type=str, required=True)
    parser.add_argument("--exp_name", type=str, required=True)
    parser.add_argument("--checkpoint", type=str, required=True)
    parser.add_argument("--lead_time", type=int, default=720)
    parser.add_argument("--n_members", type=int, required=True)
    parser.add_argument("--savetimes", nargs="+", type=int, default=[0, 6, 12, 18])
    parser.add_argument("--variables", nargs="+", default=["2t", "10u", "10v", "z_500"])
    parser.add_argument("--mask", type=str, default=None)
    args = parser.parse_args()

    if args.mask is not None:
        with open(args.mask, "rb") as f:
            mask = pickle.load(f)
    else:
        mask = None

    run_aifs(
        args.input_path,
        args.output_path,
        args.date_str,
        args.exp_name,
        args.checkpoint,
        args.lead_time,
        args.variables,
        args.savetimes,
        args.n_members,
        mask
    )
