import customtkinter as ctk
import simulator

fixed_update_interval, slow_fixed_update_interval = 0.01, 1


def gui():
    # styling
    ctk.set_appearance_mode("dark")
    ctk.set_default_color_theme("green")
    
    # window
    global window
    window = ctk.CTk()
    window.geometry("800x600")
    window.resizable(False, False)
    
    # construct frames
    sidebar_frame = ctk.CTkFrame(master=window, height=600, width=270, corner_radius=0)
    sidebar_frame.place(x=0, y=0)
    main_frame = ctk.CTkFrame(master=window, height=600, width=530, corner_radius=0, fg_color="#171717")
    main_frame.place(x=270, y=0)
    
    
    
    logo_label = ctk.CTkLabel(master=sidebar_frame, text="MOLECULAR\nSIMULATOR", font=ctk.CTkFont(family="Avant Garde", size=30, weight="bold"))
    logo_label.place(x=40, y=26)
    
    run_button = ctk.CTkButton(master=sidebar_frame, text="RUN SIMULATION", width=200, height=42, font=ctk.CTkFont(family="Tahoma", size=20, weight="bold"), command=run_simulation_callback)
    run_button.place(x=35, y=120)
    
    disclaimer_label = ctk.CTkLabel(master=sidebar_frame, text="DISCLAIMER:\nsometimes the gfx backend\nmay crash and the application\nwill need to be restarted", font=ctk.CTkFont(family="Tahoma", size=15))
    disclaimer_label.place(x=20, y=480)
    
    molecule_amount_slider = ctk.CTkSlider(master=main_frame, width=300, height=20, number_of_steps=149, command=molecule_amount_slider_callback, from_=1, to=150)
    molecule_amount_slider.place(x=20, y=20)
    molecule_amount_slider.set(50)
    simulator.molecule_amount = 50
    global molecule_amount_label
    molecule_amount_label = ctk.CTkLabel(master=main_frame, width=200, height=20, text="MOLECULE AMOUNT [50]", font=ctk.CTkFont(family="Tahoma", weight="bold"))
    molecule_amount_label.place(x=320, y=20)
    
    fixed_update_interval_slider = ctk.CTkSlider(master=main_frame, width=300, height=20, number_of_steps=60, command=fixed_update_interval_slider_callback, from_=0.01, to=10)
    fixed_update_interval_slider.place(x=20, y=50)
    fixed_update_interval_slider.set(1)
    global fixed_update_interval_label
    fixed_update_interval_label = ctk.CTkLabel(master=main_frame, width=200, height=20, text="UPDATE INTERVAL [1 fs]", font=ctk.CTkFont(family="Tahoma", weight="bold"))
    fixed_update_interval_label.place(x=320, y=50)
    
    global randomize_switch
    randomize_switch = ctk.CTkSwitch(master=main_frame, width=300, height=20, text="", switch_height=20, switch_width=80)
    randomize_switch.place(x=240, y=80)
    randomize_switch.select(1)
    global randomize_switch_label
    randomize_switch_label = ctk.CTkLabel(master=main_frame, width=200, height=20, text="RANDOMIZE POSITIONS", font=ctk.CTkFont(family="Tahoma", weight="bold"))
    randomize_switch_label.place(x=320, y=80)

    simulation_scale_slider = ctk.CTkSlider(master=main_frame, width=300, height=20, number_of_steps=30, command=simulation_scale_slider_callback, from_=1, to=30)
    simulation_scale_slider.place(x=20, y=110)
    simulation_scale_slider.set(15)
    global simulation_scale_label
    simulation_scale_label = ctk.CTkLabel(master=main_frame, width=200, height=20, text="SIMULATION SCALE [15]", font=ctk.CTkFont(family="Tahoma", weight="bold"))
    simulation_scale_label.place(x=320, y=110)
    
    simulation_radius_slider = ctk.CTkSlider(master=main_frame, width=300, height=20, number_of_steps=17, command=simulation_radius_slider_callback, from_=4, to=20)
    simulation_radius_slider.place(x=20, y=140)
    simulation_radius_slider.set(17)
    global simulation_radius_label
    simulation_radius_label = ctk.CTkLabel(master=main_frame, width=200, height=20, text="SIMULATION RADIUS [289]", font=ctk.CTkFont(family="Tahoma", weight="bold"))
    simulation_radius_label.place(x=320, y=140)
    
    global fixed_deltatime_switch
    fixed_deltatime_switch = ctk.CTkSwitch(master=main_frame, width=300, height=20, text="", switch_height=20, switch_width=80)
    fixed_deltatime_switch.place(x=240, y=170)
    fixed_deltatime_switch.select(1)
    global fixed_deltatime_label
    fixed_deltatime_label = ctk.CTkLabel(master=main_frame, width=200, height=20, text="FIXED DELTATIME", font=ctk.CTkFont(family="Tahoma", weight="bold"))
    fixed_deltatime_label.place(x=320, y=170)
    
    # run app
    window.mainloop()
    
def run_simulation_callback():
    global fixed_update_interval, slow_fixed_update_interval, randomize_switch, fixed_deltatime_switch
    simulator.main(interval=fixed_update_interval, slow_interval=slow_fixed_update_interval, randomize=randomize_switch.get(), fixed_deltatime=fixed_deltatime_switch.get())
    
def simulation_scale_slider_callback(value):
    simulator.scale_multiplier = round(value)
    global simulation_scale_label
    simulation_scale_label.configure(text="SIMULATION SCALE ["+str(round(31-value))+"]")
    
def molecule_amount_slider_callback(value):
    simulator.molecule_amount = round(value)
    global molecule_amount_label
    molecule_amount_label.configure(text="MOLECULE AMOUNT ["+str(round(value))+"]")
    
def fixed_update_interval_slider_callback(value):
    global fixed_update_interval_label, fixed_update_interval
    fixed_update_interval = round(value**2, 12)*1e-15
    fixed_update_interval_label.configure(text="UPDATE INTERVAL ["+str(round(value**2, 12))+" fs]")
    
def simulation_radius_slider_callback(value):
    global simulation_radius_label
    simulator.molecule_sim_radius = round(value**2)
    simulation_radius_label.configure(text="SIMULATION RADIUS ["+str(round(value**2))+"]")
    
    
if __name__ == '__main__':
    gui()