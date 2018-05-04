import csv,sys
import os, datetime
import time
import numpy, math
from decimal import Decimal
from math import pi, cos, sin
from datetime import datetime, time, timedelta
cities=["Ahmedabad","Allahabad","Amritsar","Aurangabad","Bangalore","Barmer","Belgaum","Bhagalpur","Bhopal","Bhubaneshwar","Bikaner","Chennai","Chirtadurg","Dehradun","Dibrugarh","Goa","Gorakhpur","Guwahati","Gwalior","Hissar","Hyderabad","Imphal","Indore","Jabalpur","Jagdelpur","Jaipur","Jaisalmer","Jamnagar","Jodhpur","Jorhat","Kolkata","Kota","Kurnool","Lucknow","Mangalore","Mumbai","Nagpur","Nellore","New Delhi","Panjim","Patna","Pune","Raipur","Rajkot","Ramagundam","Ranchi","Ratnagiri","Raxaul","Saharanpur","Shillong","Sholapur","Sundernagar","Surat","Tezpur","Tiruchirapalli","Trivandram","Veraval","Vishakhapatnam"]
DBT=[]
dew_point=[]
RH=[]
atm_prs=[]
global_horizontal_radiation=[]
direct_normal=[]
diffuse_horizontal=[]
wind_speed=[]
day_number=[]
B=[]
declination=[]
hour_angle=[]
date_lst=[]
time_lst=[]
avg_dbt=[]
occu_heat_gain=[]
equip_heat_gain=[]
light_heat_gain=[]
solar_time_minute=[]
solar_time_in_minute=[]
dbt_seven_day_mean=[]
dbt_thirty_day_mean=[]
daily_mean=[]
zenith_angle=[]
incident_angle_south=[]
incident_angle_west=[]
incident_angle_east=[]
incident_angle_north=[]
incident_angle_horizontal=[]
solar_radiation_south=[]
solar_radiation_west=[]
solar_radiation_east=[]
solar_radiation_north=[]
solar_radiation_horizontal=[]
dbt_temp=[]
ASHRAE_std_55=[]
ASHRAE_Adaptive=[]
IMAC_MM_NT=[]
IMAC_MM_UL=[]
IMAC_MM_LL=[]
set_point=[]
set_point_temperature=[]
radiant_gain_south_BAU=[]
radiant_gain_south_ECBC=[]
radiant_gain_west_BAU=[]
radiant_gain_west_ECBC=[]
radiant_gain_east_BAU=[]
radiant_gain_east_ECBC=[]
radiant_gain_north_BAU=[]
radiant_gain_north_ECBC=[]
radiant_gain_horizontal_BAU=[]
radiant_gain_horizontal_ECBC=[]
total_internal_gain_bau=[]
total_internal_gain_ecbc=[]
sol_air_temp_south=[]
sol_air_temp_west=[]
sol_air_temp_east=[]
sol_air_temp_north=[]
sol_air_temp_roof=[]
heat_gain_south_bau=[]
heat_gain_south_ecbc=[]
heat_gain_west_bau=[]
heat_gain_west_ecbc=[]
heat_gain_east_bau=[]
heat_gain_east_ecbc=[]
heat_gain_north_bau=[]
heat_gain_north_ecbc=[]
heat_gain_roof_bau=[]
heat_gain_roof_ecbc=[]
envelope_gain_bau=[]
envelope_gain_ecbc=[]
Out_humdrat=[]
In_humdrat=[]
time=[]
avg_night_time=[]
avg_night_daily=[]
occupancy_schedule=[]
lighting_schedule=[]
equipment_schedule=[]
mitigation_ECBC=[]
occupant_moisture_gain=[]
supply_air_temp_BAU=[]
supply_air_temp_ECBC=[]
supply_humdrat_bau=[]
supply_humdrat_ecbc=[]
on_coil_moisture_bau=[]
on_coil_moisture_ecbc=[]
Room_air_dew_point_bau=[]
Room_air_dew_point_ecbc=[]
infiltration=[]
envelope_temperature_rise_BAU=[]
envelope_temperature_rise_ECBC=[]
mitigation_at_night_BAU=[]
mitigation_at_night_ECBC=[]
latent_temp_rise_bau=[]
latent_temp_rise_ecbc=[]
internal_temp_rise_bau=[]
internal_temp_rise_ecbc=[]
bpt_bau=[]
bpt_ecbc=[]
CDH_bau=[]
CDH_ECBC=[]
RCDH_comfcool_bau=[]
RCDH_comfcool_ecbc=[]
WBT=[]
Tw=[]
RCDH_precool_bau=[]
RCDH_precool_ecbc=[]
RCDH_switch_bau=[]
RCDH_switch_ecbc=[]
avg_month_single=[]
avg_monthly_value=[]
month={'Jan':31, 'Feb':28, 'Mar':31, 'Apr':30, 'May':31, 'Jun':30, 'Jul':31, 'Aug':31, 'Sept':30, 'Oct':31, 'Nov':30, 'Dec':31}
abc=[]
avg_seasonal_value=[]
seas_avg=[]
min_jan=[]
min_feb=[]
min_mar=[]
min_apr=[]
min_may=[]
min_jun=[]
min_jul=[]
min_aug=[]
min_sep=[]
min_oct=[]
min_nov=[]
min_dec=[]
month_min=[]
monthly_min=[]
min_seasonal=[]
minimum_seasonal=[]
daily_max=[]
daily_min=[]
TinmaxECBC=[]
TinminECBC=[]
TinmaxBAU=[]
TinminBAU=[]
average_DBT=[]
daily_maximum=[]
daily_minimum=[]
amplitude_bau=[]
amplitude_ecbc=[]
sine_theeta=[]
angle_radian=[]
Tin_hourly_bau=[]
Tin_hourly_ecbc=[]
RCDH_night_bau=[]
RCDH_night_ecbc=[]
evap_precool_energy_hourly_bau=[]
evap_precool_energy_hourly_ecbc=[]
evap_cool_energy_hourly_bau=[]
evap_cool_energy_hourly_ecbc=[]
nightcool_energy_hourly_bau=[]
nightcool_energy_hourly_ecbc=[]
evap_temp=[]


def open_file_to_read(file_name):
    file = open(file_name, "r")
    return file
		
def open_file_to_write(file_name):
    file = open(file_name, "w+")
    return file

def find_latitude_longitude(file1,city_name):
    latitude=0
    longitude=0
    roof_U_value=0 
    wall_U_value=0
    window_U_value=0
    window_SHGC_nonnorth=0
    window_SHGC_north=0
    #print (city_name)    
    next(file1)
    for line in file1:

        row_heading=line.split(",")
        if row_heading[0]==city_name:
            latitude=row_heading[2]
            longitude=row_heading[3]
            roof_U_value=row_heading[4]
            wall_U_value=row_heading[5]
            window_U_value=row_heading[6]
            window_SHGC_nonnorth=row_heading[7]
            window_SHGC_north=row_heading[8]
            break
    return latitude, longitude, roof_U_value, wall_U_value, window_U_value, window_SHGC_nonnorth, window_SHGC_north


def building_geometry(roof_U_value, wall_U_value, building_area, orientation, building_height, num_floor, asp_ratio, WWR):
    if orientation=="E-W" or "e-w":
        Twan=int(asp_ratio)*math.sqrt(int(building_area)/int(asp_ratio))*int(building_height)
    else:
        Twan=math.sqrt(int(building_area)/int(asp_ratio))*int(building_height)
    Twae=int(Twan)/int(asp_ratio)
    Twaw=int(Twae)
    Twas=int(Twan)
    Win_area_north=float(WWR)*int(Twan)
    Win_area_east=float(WWR)*int(Twae)
    Win_area_west=float(WWR)*int(Twaw)
    Win_area_south=float(WWR)*int(Twas)
    Wall_area_north=int(Twan)-int(Win_area_north)
    Wall_area_east=int(Twae)-int(Win_area_east)
    Wall_area_west=int(Twaw)-int(Win_area_west)
    Wall_area_south=int(Twas)-int(Win_area_south)
    Roof_area=int(building_area)/int(num_floor)
    VoP_north=int(Wall_area_north)*0.015
    VoP_east=int(Wall_area_east)*0.015
    VoP_west=int(Wall_area_west)*0.015
    VoP_south=int(Wall_area_south)*0.015
    VoP_roof=int(Roof_area)*0.15
    heat_carrying_capacity_of_air=int(building_area)*0.00624 #in kW/K
    specific_heat_wall=970  #in kJ/kg K
    specific_heat_roof=480  #in kJ/kg K
    wall_density=1880       #in kg/m3
    roof_density=2427       #in kg/m3
    roof_surf_refl=0.8
    wall_surf_refl=0.6
    overall_heat_transfer_coefficient_BAU=((float(Wall_area_north+ Wall_area_west+ Wall_area_east+ Wall_area_south)*1.772)+(float(Roof_area)*2.942)+(0.3*float(building_area)*float(building_height)/7))/1000
    overall_heat_transfer_coefficient_ECBC=((float(Wall_area_north+ Wall_area_west+ Wall_area_east+ Wall_area_south)*float(wall_U_value))+(float(Roof_area)*float(roof_U_value))+(0.3*float(building_area)*float(building_height)/7))/1000
    thermal_capacitance=((float(VoP_north+VoP_east+VoP_west+VoP_south)*float(specific_heat_wall)*float(wall_density))+(float(VoP_roof)*float(specific_heat_roof)*float(roof_density)))/1000
    time_constant_BAU=thermal_capacitance/(3600*float(overall_heat_transfer_coefficient_BAU))
    time_constant_ECBC=thermal_capacitance/(3600*float(overall_heat_transfer_coefficient_ECBC))
    
    return Twan, Twae, Twaw, Twas, Win_area_north, Win_area_east, Win_area_west, Win_area_south, Wall_area_north, Wall_area_east, Wall_area_west, Wall_area_south, Roof_area, VoP_north, VoP_east, VoP_west, VoP_south, heat_carrying_capacity_of_air, specific_heat_wall, specific_heat_roof, roof_surf_refl, wall_surf_refl, time_constant_BAU, time_constant_ECBC, thermal_capacitance

def read_climate_variables(f1):
    next(f1)
    
    for line in f1:
        row_heading=line.split(",")
        DBT_value=row_heading[2]	
        DBT.append(DBT_value)
        dew_point_value=row_heading[3]
        dew_point.append(dew_point_value)
        RH_value=row_heading[4]
        RH.append(RH_value)
        atm_prs_value=row_heading[5]
        atm_prs.append(atm_prs_value)
        global_horizontal_radiation_value=row_heading[6]
        global_horizontal_radiation.append(global_horizontal_radiation_value)
        direct_normal_value=row_heading[7]
        direct_normal.append(direct_normal_value)
        diffuse_horizontal_value=row_heading[8]
        diffuse_horizontal.append(diffuse_horizontal_value)
        wind_speed_value=row_heading[9]
        wind_speed.append(wind_speed_value)

    return DBT, dew_point, RH,atm_prs,global_horizontal_radiation,direct_normal, diffuse_horizontal, wind_speed     

def daily_average(file_comfort, DBT):
    i=0
    j=0
    length_dbt=len(DBT)
    dbt_value=0
    
    while(j<length_dbt):
      
        dbt_value=dbt_value+float(DBT[i])
        if i!=0:
            if i%24==23:
                avg_dbt_value=float(dbt_value)/24
                avg_dbt.append(avg_dbt_value)
                dbt_value=0
        i=i+1
        j=j+1
        
    return avg_dbt

def avg_night(DBT):
    
    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\avg_night.txt")
    avg_night_file=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\avg_night.txt")
    
    i=0
    j=0
    length_dbt=len(DBT)
    dbt_value=0
    avg_night_value=0
    
    while(j<length_dbt):

        if i<6 or i>17:

            dbt_value=dbt_value+float(DBT[j])
            if i!=1:
                
                if i%24==23:
                    
                    avg_night_value=float(dbt_value)/12
                    avg_night_time.append(avg_night_value)
                    
                    avg_night_file.write(str(avg_night_value)+"\n")
                    dbt_value=0
                    i=-1
                    
        i=i+1
        j=j+1
    return avg_night_time
    avg_night_file.close()

    
def daily_avg_night(avg_night):

    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\avg_night_daily.txt")
    avg_night_daily_file=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\avg_night_daily.txt")
    
    avg_night_file2=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\avg_night_daily.txt")
    
    for items in avg_night_file2:
        avg_night_time.append(items)
        
    for items in avg_night_time:
   
        i=0
        while i<24:
            
            avg_night_value=items
                
            avg_night_daily_file.write(str(avg_night_value)+"\n")
            i=i+1
            
    avg_night_daily_file.close()      
    

def seven_day_mean(avg_dbt):
    i=0
    for items in avg_dbt:
        day_7=avg_dbt[i-6]
        day_6=avg_dbt[i-5]
        day_5=avg_dbt[i-4]
        day_4=avg_dbt[i-3]
        day_3=avg_dbt[i-2]
        day_2=avg_dbt[i-1]
        day_1=avg_dbt[i]
        
        dbt_seven_day_mean_value=((0.33*day_1)+(0.23*day_2)+(0.16*day_3)+(0.11*day_4)+(0.08*day_5)+(0.08*day_6)+(0.04*day_7))
        dbt_seven_day_mean.append(dbt_seven_day_mean_value)
        i=i+1
        
    return dbt_seven_day_mean

def thirty_day_mean(daily_mean_IMAC):
    
    i=0
    for items in daily_mean_IMAC:
        day_30=avg_dbt[i-29]
        day_29=avg_dbt[i-28]
        day_28=avg_dbt[i-27]
        day_27=avg_dbt[i-26]
        day_26=avg_dbt[i-25]
        day_25=avg_dbt[i-24]
        day_24=avg_dbt[i-23]
        day_23=avg_dbt[i-22]
        day_22=avg_dbt[i-21]
        day_21=avg_dbt[i-20]
        day_20=avg_dbt[i-19]
        day_19=avg_dbt[i-18]
        day_18=avg_dbt[i-17]
        day_17=avg_dbt[i-16]
        day_16=avg_dbt[i-15]
        day_15=avg_dbt[i-14]
        day_14=avg_dbt[i-13]
        day_13=avg_dbt[i-12]
        day_12=avg_dbt[i-11]
        day_11=avg_dbt[i-10]
        day_10=avg_dbt[i-9]
        day_09=avg_dbt[i-8]
        day_08=avg_dbt[i-7]
        day_07=avg_dbt[i-6]
        day_06=avg_dbt[i-5]
        day_05=avg_dbt[i-4]
        day_04=avg_dbt[i-3]
        day_03=avg_dbt[i-2]
        day_02=avg_dbt[i-1]
        day_01=avg_dbt[i]
        
        dbt_thirty_day_mean_value=(float(day_01)+float(day_02)+float(day_03)+float(day_04)+float(day_05)+float(day_06)+float(day_07)+float(day_08)+float(day_09)+float(day_10)+float(day_11)+float(day_12)+float(day_13)+float(day_14)+float(day_15)+float(day_16)+float(day_17)+float(day_18)+float(day_19)+float(day_20)+float(day_21)+float(day_22)+float(day_23)+float(day_24)+float(day_25)+float(day_26)+float(day_27)+float(day_28)+float(day_29)+float(day_30))/(30)
        dbt_thirty_day_mean.append(dbt_thirty_day_mean_value)
        
        i=i+1
        
    return dbt_thirty_day_mean

def daily_mean_IMAC(DBT):
    i=0
    j=0
    length_DBT=len(DBT)
    
    while j<length_DBT:
        
        i=0
        
        del dbt_temp[:]
        
        while i<24:
            
            dbt_temp.append(DBT[j])
            i=i+1
            j=j+1
            
        daily_max_value=max(dbt_temp)
        
        daily_min_value=min(dbt_temp)
        
        daily_mean_value=(float(daily_max_value)+float(daily_min_value))/2

        daily_max.append(daily_max_value)

        daily_min.append(daily_min_value)
        
        daily_mean.append(daily_mean_value)

    return daily_mean, daily_max, daily_min

def comfort_model(avg_dbt, dbt_seven_day_mean, dbt_thirty_day_mean, daily_mean):

    i=0
    while i<365:
        ASHRAE_std_55_value=22.5
        ASHRAE_std_55.append(ASHRAE_std_55_value)
        i=i+1
        
    for items in dbt_seven_day_mean:
        ASHRAE_Adaptive_value=(0.31*float(items))+17.8
        ASHRAE_Adaptive.append(ASHRAE_Adaptive_value)

    for items in dbt_thirty_day_mean:
        IMAC_MM_NT_value=0.28*float(items)+17.87
        IMAC_MM_NT.append(IMAC_MM_NT_value)

    for items in IMAC_MM_NT:
        IMAC_MM_UL_value=float(items)+3.46
        IMAC_MM_UL.append(IMAC_MM_UL_value)
        
        IMAC_MM_LL_value=float(items)-3.46
        IMAC_MM_LL.append(IMAC_MM_LL_value)
        
    return ASHRAE_std_55, ASHRAE_Adaptive, IMAC_MM_NT, IMAC_MM_UL, IMAC_MM_LL

def set_point_temp(comf_temp_model, DBT):

    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\comfort_model.txt")
    comfort_model_file=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\comfort_model.txt")
    
    file_internal=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\Internal_load.csv")

    next(file_internal)
    
    for line in file_internal:
        row_heading=line.split(",")
        occupancy_schedule.append(row_heading[3])

    for items in comf_temp_model:
        
        i=0
        while i<24:
            if float(occupancy_schedule[i])>0:
                set_point_value=items
            
                
            else:
                set_point_value=DBT[i]
            set_point_value=round(float(set_point_value),2)
            set_point.append(set_point_value)
            comfort_model_file.write(str(set_point_value)+"\n")
            i=i+1
    comfort_model_file.close()
    file_internal.close()
    return set_point
        
def part_press(atm_prs,W):
    i=0
    for items in atm_prs:
        Part_press_value = atm_prs[i]*W/(0.62198+W)
        Part_press.append(Part_press_value)
        i+1
    return Part_press

def humdity_ratio(DBT, RH, atm_prs, set_point):

    comfort_model_file=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\comfort_model.txt")
    
        
    C1 = -5674.5359
    
    C2 = 6.3925247
    
    C3 = -0.009677843
    
    C4 = 0.00000062215701
    
    C5 = 2.0747825E-09
    
    C6 = -9.484024E-13
    
    C7 = 4.1635019
    
    C8 = -5800.2206
    
    C9 = 1.3914993
    
    C10 = -0.048640239
    
    C11 = 0.000041764768
    
    C12 = -0.000000014452093
    
    C13 = 6.5459673

    i=0
    for items in DBT:
        
        Tk=float(DBT[i])+273.5
        
        if Tk<=273.5:
            
            Pws_value = math.exp(C1/Tk + C2 + C3*Tk + C4*Tk**2 + C5*Tk**3 + C6*Tk**4 + C7*math.log(Tk)) / 1000

        else:

            Pws_value = math.exp(C8/Tk + C9 + C10*Tk + C11*Tk**2 + C12*Tk**3 + C13*math.log(Tk)) / 1000
            
        Out_humdrat_value=0.621945*(float(RH[i])/100)*float(Pws_value)/((float(atm_prs[i])/1000)-((float(RH[i])/100)*float(Pws_value)))
        Out_humdrat_value = round(Out_humdrat_value,4)
        Out_humdrat.append(Out_humdrat_value)
        #print("our_humdrat"+str(Out_humdrat_value)+"\n")
        i=i+1
        
    i=0
    for items in set_point:
        Tki=float(set_point[i])+273.5
        
        if Tki<=273.5:
            
            Pwsi_value = math.exp(C1/Tki + C2 + C3*Tki+ C4*Tki**2 + C5*Tki**3 + C6*Tki**4 + C7*math.log(Tki)) / 1000

        else:

            Pwsi_value = math.exp(C8/Tki + C9 + C10*Tki + C11*Tki**2 + C12*Tki**3 + C13*math.log(Tki)) / 1000
            
        In_humdrat_value=0.621945*0.6*float(Pwsi_value)/((float(atm_prs[i])/1000)-0.6*float(Pwsi_value))
        In_humdrat_value = round(In_humdrat_value,4)
        In_humdrat.append(In_humdrat_value)
        #print("in_humdrat"+str(In_humdrat_value)+"\n")
        i=i+1
    
    return Out_humdrat, In_humdrat
    comfort_model_file.close()
   
def Hum_rat(Tdb, RH, P):

#Function to calculate humidity ratio [kg H2O/kg air]
#Given dry bulb and wet bulb temperature inputs [degC]
#ASHRAE Fundamentals handbood (2005)
#Tdb = Dry bulb temperature [degC]
#RH = Relative Humidity [Fraction or %]
#P = Ambient Pressure [kPa]

    Pws = Sat_press(Tdb)

    result = 0.62198*0.6*Pws/(P - RH*Pws)    # Equation 22, 24, p6.8

    return result

def wet_bulb(DBT, RH, atm_prs):
    
    W_normal = Hum_rat2(DBT, RH, atm_prs)

    result = Tdb

    #Solves to within 0.001% accuracy using Newton-Rhapson   

    W_new = Hum_rat(DBT, result, atm_prs)

    i=0

    for items in DBT:
        

        while abs((W_new - W_normal) / W_normal) > 0.00001:

            W_new2 = Hum_rat(DBT, result - 0.001, atm_prs)

            dw_dtwb = (W_new - W_new2) / 0.001

            result = result - (W_new - W_normal) / dw_dtwb

            W_new = Hum_rat(DBT, result, atm_prs)

    return result

def wet_bulb(DBT,RH):
    
    i=0
    for items in DBT:
        
        
        Tw_value=(float(DBT[i])*(math.atan(0.151977*((float(RH[i]))**0.5))))+(math.atan(float(DBT[i])+float(RH[i])))-(math.atan(float(RH[i])-1.676331))+((0.00391838*(float(RH[i])**1.5))*math.atan(0.023*float(RH[i])))-4.686035
        Tw.append(Tw_value)
        i=i+1
    return Tw

def Dew_point(atm_prs, W):

    #Function to compute the dew point temperature (deg C)

    #From page 6.9 equation 39 and 40 in ASHRAE Fundamentals handbook (2005)

    #P = ambient pressure [kPa]

    #W = humidity ratio [kg/kg dry air]

    #Valid for Dew Points less than 93 C

    C14 = 6.54

    C15 = 14.526

    C16 = 0.7389

    C17 = 0.09486

    C18 = 0.4569

    i=0
    for items in atm_prs:

        Pw = Part_press(atm_prs, W)

        alpha = math.log(Pw)

        Tdp1 = C14 + C15*alpha + C16*alpha**2 + C17*alpha**3 + C18*Pw**0.1984

        Tdp2 = 6.09 + 12.608*alpha + 0.4959*alpha**2

        if Tdp1 >= 0:

            result = Tdp1

        else:

            result = Tdp2
            
        i=i+1

    return result
   
def convert_local_Solar(longitude, file_Fabric):
    next(file_Fabric)
    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\solar_time_in_minute.txt")
    solar_time_in_minute_file=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\solar_time_in_minute.txt")
    for line in file_Fabric:
        row_heading=line.split(",")

        date_value=row_heading[0]
        date_lst.append(date_value)
        time_value=row_heading[1]+":00"
        time_lst.append(time_value)
        
    longitude=float(longitude)
    i=0
    solar_time=0
    
    for items in date_lst:
		
        year=date_lst[i].split("-")[2]
        month=date_lst[i].split("-")[1]
        day=date_lst[i].split("-")[0]
        date_in_ddmmyyyy=day+"-"+month+"-"+year
        
		
        hour=int(time_lst[i].split(":")[0])
        minute=time_lst[i].split(":")[1]
        sec=time_lst[i].split(":")[2]
        if hour<24:
            hour=time_lst[i].split(":")[0]
            time_in_hhmmss=str(hour)+":"+str(minute)+":"+str(sec)
            
            i=i+1
            dt = datetime.strptime(date_in_ddmmyyyy + ' ' + time_in_hhmmss, '%d-%m-%Y %H:%M:%S')	
            
            gamma = 2 * pi / 365 * (dt.timetuple().tm_yday + 1 + float(dt.hour - 12) / 24)
		
            eqtime = 229.18 * (0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) - 0.014615 * cos(2 * gamma) - 0.040849 * sin(2 * gamma))
		
            decl = 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma)- 0.006758 * cos(2 * gamma) + 0.000907 * sin(2 * gamma) - 0.002697 * cos(3 * gamma) + 0.00148 * sin(3 * gamma)
    		
            time_offset = eqtime + 4 *(82.5-longitude)
            solar_time_in_minute_value = dt.hour * 60 + dt.minute + dt.second / 60 + time_offset
            solar_time_in_minute_file.write(str(solar_time_in_minute_value)+"\n")

    solar_time_in_minute_file.close()


def monthly_average(DBT):
    
    DBT_jan=0
    i=0
    while i<745:
        DBT_jan=DBT_jan+float(DBT[i])
        i=i+1
    avg_jan = DBT_jan/744

    DBT_feb=0
    i=745
    while i<1417:
        DBT_feb=DBT_feb+float(DBT[i])
        i=i+1
    avg_feb = DBT_feb/672
    
    DBT_mar=0
    i=1417
    while i<2161:
        DBT_mar=DBT_mar+float(DBT[i])
        i=i+1
    avg_mar = DBT_mar/744

    DBT_apr=0
    i=2161
    while i<2881:
        DBT_apr=DBT_apr+float(DBT[i])
        i=i+1
    avg_apr = DBT_apr/720

    DBT_may=0
    i=2881
    while i<3625:
        DBT_may=DBT_may+float(DBT[i])
        i=i+1
    avg_may = DBT_may/744

    DBT_jun=0
    i=3625
    while i<4345:
        DBT_jun=DBT_jun+float(DBT[i])
        i=i+1
    avg_jun = DBT_jun/720

    DBT_jul=0
    i=4345
    while i<5089:
        DBT_jul=DBT_jul+float(DBT[i])
        i=i+1
    avg_jul = DBT_jul/744

    DBT_aug=0
    i=5089
    while i<5833:
        DBT_aug=DBT_aug+float(DBT[i])
        i=i+1
    avg_aug = DBT_aug/744

    DBT_sep=0
    i=5833
    while i<6533:
        DBT_sep=DBT_sep+float(DBT[i])
        i=i+1
    avg_sep = DBT_sep/720

    DBT_oct=0
    i=6533
    while i<7273:
        DBT_oct=DBT_oct+float(DBT[i])
        i=i+1
    avg_oct = DBT_oct/744

    DBT_nov=0
    i=7273
    while i<8017:
        DBT_nov=DBT_nov+float(DBT[i])
        i=i+1
    avg_nov = DBT_nov/720

    DBT_dec=0
    i=8017
    while i<=8759:
        DBT_dec=DBT_dec+float(DBT[i])
        i=i+1
    avg_dec = DBT_dec/744
    
    avg_month_single.append(avg_jan)
    avg_month_single.append(avg_feb)
    avg_month_single.append(avg_mar)
    avg_month_single.append(avg_apr)
    avg_month_single.append(avg_may)
    avg_month_single.append(avg_jun)
    avg_month_single.append(avg_jul)
    avg_month_single.append(avg_aug)
    avg_month_single.append(avg_sep)
    avg_month_single.append(avg_oct)
    avg_month_single.append(avg_nov)
    avg_month_single.append(avg_dec)

    return avg_month_single, avg_jan, avg_feb, avg_mar, avg_apr, avg_may, avg_jun, avg_jul, avg_aug, avg_sep, avg_oct, avg_nov, avg_dec

def avg_monthly (avg_month_single):
    i=0
    while i<month["Jan"]:
        avg_monthly_value.append(float(avg_month_single[0]))
        i=i+1

    i=0
    while i<month["Feb"]:
        avg_monthly_value.append(float(avg_month_single[1]))
        i=i+1

    i=0
    while i<month["Mar"]:
        avg_monthly_value.append(float(avg_month_single[2]))
        i=i+1

    i=0
    while i<month["Apr"]:
        avg_monthly_value.append(float(avg_month_single[3]))
        i=i+1

    i=0
    while i<month["May"]:
        avg_monthly_value.append(float(avg_month_single[4]))
        i=i+1

    i=0
    while i<month["Jun"]:
        avg_monthly_value.append(float(avg_month_single[5]))
        i=i+1

    i=0
    while i<month["Jul"]:
        avg_monthly_value.append(float(avg_month_single[6]))
        i=i+1

    i=0
    while i<month["Aug"]:
        avg_monthly_value.append(float(avg_month_single[7]))
        i=i+1

    i=0
    while i<month["Sept"]:
        avg_monthly_value.append(float(avg_month_single[8]))
        i=i+1

    i=0
    while i<month["Oct"]:
        avg_monthly_value.append(float(avg_month_single[9]))
        i=i+1

    i=0
    while i<month["Nov"]:
        avg_monthly_value.append(float(avg_month_single[10]))
        i=i+1

    i=0
    while i<month["Dec"]:
        avg_monthly_value.append(float(avg_month_single[11]))
        i=i+1

    return avg_monthly_value

def monthly_average_temparature(avg_monthly_value):
    j=0
    
    while j<365:
        i=0
        while i<24:
            abc_value=float(avg_monthly_value[j])
            abc.append(avg_monthly_value[j])
            i=i+1
        j=j+1
        
    return abc
        
def seasonal_average(avg_jan, avg_feb, avg_mar, avg_apr, avg_may, avg_jun, avg_jul, avg_aug, avg_sep, avg_oct, avg_nov, avg_dec, DBT):
    avg_season_1 = (float(avg_dec)+float(avg_jan)+float(avg_feb))/3
    avg_season_2 = (float(avg_mar)+float(avg_apr)+float(avg_may))/3
    avg_season_3 = (float(avg_jun)+float(avg_jul)+float(avg_aug))/3
    avg_season_4 = (float(avg_sep)+float(avg_oct)+float(avg_nov))/3
    
    i=0
    while i<month["Jan"]:
        avg_seasonal_value.append(float(avg_season_1))
        i=i+1

    i=0
    while i<month["Feb"]:
        avg_seasonal_value.append(float(avg_season_1))
        i=i+1

    i=0
    while i<month["Mar"]:
        avg_seasonal_value.append(float(avg_season_2))
        i=i+1

    i=0
    while i<month["Apr"]:
        avg_seasonal_value.append(float(avg_season_2))
        i=i+1

    i=0
    while i<month["May"]:
        avg_seasonal_value.append(float(avg_season_2))
        i=i+1

    i=0
    while i<month["Jun"]:
        avg_seasonal_value.append(float(avg_season_3))
        i=i+1

    i=0
    while i<month["Jul"]:
        avg_seasonal_value.append(float(avg_season_3))
        i=i+1

    i=0
    while i<month["Aug"]:
        avg_seasonal_value.append(float(avg_season_3))
        i=i+1

    i=0
    while i<month["Sept"]:
        avg_seasonal_value.append(float(avg_season_4))
        i=i+1

    i=0
    while i<month["Oct"]:
        avg_seasonal_value.append(float(avg_season_4))
        i=i+1

    i=0
    while i<month["Nov"]:
        avg_seasonal_value.append(float(avg_season_4))
        i=i+1

    i=0
    while i<month["Dec"]:
        avg_seasonal_value.append(float(avg_season_1))
        i=i+1

    j=0
    
    while j<365:
        i=0
        while i<24:
            seas_avg_value=float(avg_seasonal_value[j])
            seas_avg.append(seas_avg_value)
            i=i+1
        j=j+1
    return seas_avg

def monthly_minimum(DBT):
    i=0
    while i<745:
        min_jan.append(float(DBT[i]))
        i=i+1
    Min_jan = min(min_jan)

    i=745
    while i<1417:
        min_feb.append(float(DBT[i]))
        i=i+1
    Min_feb = min(min_feb)
    
    i=1417
    while i<2161:
        min_mar.append(float(DBT[i]))
        i=i+1
    Min_mar = min(min_mar)

    i=2161
    while i<2881:
        min_apr.append(float(DBT[i]))
        i=i+1
    Min_apr = min(min_apr)

    i=2881
    while i<3625:
        min_may.append(float(DBT[i]))
        i=i+1
    Min_may = min(min_may)

    i=3625
    while i<4345:
        min_jun.append(float(DBT[i]))
        i=i+1
    Min_jun = min(min_jun)

    i=4345
    while i<5089:
        min_jul.append(float(DBT[i]))
        i=i+1
    Min_jul = min(min_jul)

    i=5089
    while i<5833:
        min_aug.append(float(DBT[i]))
        i=i+1
    Min_aug = min(min_aug)

    i=5833
    while i<6533:
        min_sep.append(float(DBT[i]))
        i=i+1
    Min_sep = min(min_sep)

    i=6533
    while i<7273:
        min_oct.append(float(DBT[i]))
        i=i+1
    Min_oct = min(min_oct)

    i=7273
    while i<8017:
        min_nov.append(float(DBT[i]))
        i=i+1
    Min_nov = min(min_nov)

    i=8017
    while i<=8759:
        min_dec.append(float(DBT[i]))
        i=i+1
    Min_dec = min(min_dec)

    i=0
    while i<month["Jan"]:
        month_min.append(float(Min_jan))
        i=i+1

    i=0
    while i<month["Feb"]:
        month_min.append(float(Min_feb))
        i=i+1

    i=0
    while i<month["Mar"]:
        month_min.append(float(Min_mar))
        i=i+1

    i=0
    while i<month["Apr"]:
        month_min.append(float(Min_apr))
        i=i+1

    i=0
    while i<month["May"]:
        month_min.append(float(Min_may))
        i=i+1

    i=0
    while i<month["Jun"]:
        month_min.append(float(Min_jun))
        i=i+1

    i=0
    while i<month["Jul"]:
        month_min.append(float(Min_jul))
        i=i+1

    i=0
    while i<month["Aug"]:
        month_min.append(float(Min_aug))
        i=i+1

    i=0
    while i<month["Sept"]:
        month_min.append(float(Min_sep))
        i=i+1

    i=0
    while i<month["Oct"]:
        month_min.append(float(Min_oct))
        i=i+1

    i=0
    while i<month["Nov"]:
        month_min.append(float(Min_nov))
        i=i+1

    i=0
    while i<month["Dec"]:
        month_min.append(float(Min_dec))
        i=i+1

    j=0
    
    while j<365:
        i=0
        while i<24:
            monthly_min_value=float(month_min[j])
            monthly_min.append(monthly_min_value)
            i=i+1
        j=j+1
        
    return Min_jan, Min_feb, Min_mar, Min_apr, Min_may, Min_jun, Min_jul, Min_aug, Min_sep, Min_oct, Min_nov, Min_dec, monthly_min

def seasonal_minimum(Min_jan, Min_feb, Min_mar, Min_apr, Min_may, Min_jun, Min_jul, Min_aug, Min_sep, Min_oct, Min_nov, Min_dec):

    season1 = min(Min_dec, Min_jan, Min_feb)
    season2 = min(Min_mar, Min_apr, Min_may)
    season3 = min(Min_jun, Min_jul, Min_aug)
    season4 = min(Min_sep, Min_oct, Min_nov)

    i=0
    while i<month["Jan"]:
        min_seasonal.append(float(season1))
        i=i+1

    i=0
    while i<month["Feb"]:
        min_seasonal.append(float(season1))
        i=i+1

    i=0
    while i<month["Mar"]:
        min_seasonal.append(float(season2))
        i=i+1

    i=0
    while i<month["Apr"]:
        min_seasonal.append(float(season2))
        i=i+1

    i=0
    while i<month["May"]:
        min_seasonal.append(float(season2))
        i=i+1

    i=0
    while i<month["Jun"]:
        min_seasonal.append(float(season3))
        i=i+1

    i=0
    while i<month["Jul"]:
        min_seasonal.append(float(season3))
        i=i+1

    i=0
    while i<month["Aug"]:
        min_seasonal.append(float(season3))
        i=i+1

    i=0
    while i<month["Sept"]:
        min_seasonal.append(float(season4))
        i=i+1

    i=0
    while i<month["Oct"]:
        min_seasonal.append(float(season4))
        i=i+1

    i=0
    while i<month["Nov"]:
        min_seasonal.append(float(season4))
        i=i+1

    i=0
    while i<month["Dec"]:
        min_seasonal.append(float(season1))
        i=i+1

    j=0
    
    while j<365:
        i=0
        while i<24:
            minimum_seasonal_value=float(min_seasonal[j])
            minimum_seasonal.append(minimum_seasonal_value)
            i=i+1
        j=j+1

    return minimum_seasonal
            

def hour_angle_calculation(file_solar_time_in_minute):
    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\hour_angle.txt")
    hour_angle_file=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\hour_angle.txt")
    hour_angle_file.write("hour_angle\n")
    for items in file_solar_time_in_minute:
        solar_time_in_minute.append(items)
    for items in solar_time_in_minute:
        if float(items)/4<0:
            hour_angle_value=float(items)/4+180
        else:
            hour_angle_value=float(items)/4-180
        hour_angle_file.write(str(hour_angle_value)+"\n")
       
    hour_angle_file.close()
    
        

def do_calculations_for_solar_radiation(solar_time_in_minute,file_Solar, hour_angle_file, latitude):
    
    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\solar_angles.txt")
    solar_angles_file=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\solar_angles.txt")
    solar_angles_file.write("zenith_angle,altitude_angle,incident_angle_south_,incident_angle_west,incident_angle_east,incident_angle_north,incident_angle_roof\n")
    
    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\B.txt")
    B_file=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\B.txt")
    B_file.write("B\n")
    
    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\declination.txt")
    declination_file=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\declination.txt")
    declination_file.write("Declination\n")
    
    next(file_Solar)
    next(hour_angle_file)
    
    for items in hour_angle_file:
        
        hour_angle.append(items)
        
    for line in file_Solar:
        
        row_heading=line.split(",")
        day_number.append(row_heading[4])
    
    i=0
    for items in day_number:
        
        B_individual_values=float(int(day_number[i])-1)*(float(360)/float(365))
        
        B.append(B_individual_values)
        B_file.write(str(B_individual_values)+"\n")
        
        i=i+1
        
    B_file.close()
    
    i=0
    for items in B:
        cosB= math.cos(math.degrees(B[i]))
        sinB=math.sin(math.degrees(B[i]))
        cos2B=math.cos(2*math.degrees(B[i]))
        sin2B=math.sin(2*math.degrees(B[i]))
        cos3B=math.cos(3*math.degrees(B[i]))
        sin3B=math.sin(3*math.degrees(B[i]))
        
        declination_individual_values=math.degrees(0.006918-(0.399912*cosB)+(0.070257*sinB)-(0.006758*cos2B)+(0.000907*sin2B)-(0.002697*cos3B)+(0.00148*sin3B))

        declination.append(declination_individual_values)
        
        declination_file.write(str(declination_individual_values)+"\n")
        
        i=i+1
    declination_file.close

    declination_file=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\declination.txt")
    next(declination_file)
    for items in declination_file:
        declination.append(items)
    i=0
    for items in hour_angle:
        zenith_angle_value1=math.sin(float(latitude))*math.sin(float(declination[i]))
        zenith_angle_value2=math.cos(float(latitude))*math.cos(float(declination[i]))*math.cos(float(items))
        zenith_angle_value3=zenith_angle_value1+zenith_angle_value2
        zenith_angle_value=math.degrees(math.acos(zenith_angle_value3))
        zenith_angle_value=round(zenith_angle_value,2)
        solar_angles_file.write(str(zenith_angle_value)+",")
        
        altitude_angle_value=float(zenith_angle_value)-90
        altitude_angle_value=round(altitude_angle_value,2)
        solar_angles_file.write(str(altitude_angle_value)+",")
        
        sin_decl=math.sin(math.radians(float(declination[i])))
        cos_decl=math.cos(math.radians(float(declination[i])))
        sin_lati=math.sin(math.radians(float(latitude)))
        cos_lati=math.cos(math.radians(float(latitude)))
        sin_hour=math.sin(math.radians(float(items)))
        cos_hour=math.cos(math.radians(float(items)))
        
        value1=sin_decl*sin_lati
        value2=sin_decl*cos_lati   
        value3=cos_decl*cos_lati*cos_hour    
        value4=cos_decl*sin_lati*cos_hour     
        value5=cos_decl*sin_hour
        
        if float(hour_angle[i])>-90 and float(hour_angle[i])<=90 :
            
            a=(float(value1)*math.cos(math.radians(90)))-(float(value2)*math.sin(math.radians(90))*math.cos(0))+(float(value3)*math.cos(math.radians(90)))+(float(value4)*math.sin(math.radians(90)))*math.cos(0)+(float(value5)*math.sin(math.radians(90)))*math.sin(0)
            acos_a=math.degrees(math.acos(a))
            acos_a=round(acos_a,0)
            
                
            b=(float(value1)*math.cos(math.radians(90)))-(float(value2)*math.sin(math.radians(90))*math.cos(math.radians(90)))+(float(value3)*math.cos(math.radians(90))*math.cos(math.radians(90)))+(float(value4)*math.sin(math.radians(90))*math.cos(math.radians(90)))+(float(value5)*math.sin(math.radians(90))*math.sin(math.radians(90)))
            acos_b=math.degrees(math.acos(b))
            acos_b=round(acos_b,0)
            
        
            c=(float(value1)*math.cos(math.radians(90)))-(float(value2)*math.sin(math.radians(90))*math.cos(math.radians(270)))+(float(value3)*math.cos(math.radians(270))*math.cos(math.radians(90)))+(float(value4)*math.sin(math.radians(90))*math.cos(math.radians(270)))+(float(value5)*math.sin(math.radians(90))*math.sin(math.radians(270)))
            acos_c=math.degrees(math.acos(c))
            acos_c=round(acos_c,0)
            
        
            d=(float(value1)*math.cos(math.radians(90)))-(float(value2)*math.sin(math.radians(90))*math.cos(math.radians(180)))+(float(value3)*math.cos(math.radians(180))*math.cos(math.radians(90)))+(float(value4)*math.sin(math.radians(90))*math.cos(math.radians(180)))+(float(value5)*math.sin(math.radians(90))*math.sin(math.radians(180)))
            acos_d=math.degrees(math.acos(d))
            acos_d=round(acos_d,0)

            e=(float(value1)*math.cos(math.radians(0)))-(float(value2)*math.sin(math.radians(0))*math.cos(math.radians(0)))+(float(value3)*math.cos(math.radians(0)))+(float(value4)*math.sin(math.radians(0))*math.cos(math.radians(0)))+(float(value5)*math.sin(math.radians(0))*math.sin(math.radians(0))) 
            acos_e=math.degrees(math.acos(e))
            acos_e=round(acos_e,0)
            
        else:
            acos_a=0
            acos_b=0
            acos_c=0
            acos_d=0
            acos_e=0
            
        solar_angles_file.write(str(acos_a)+",")
        solar_angles_file.write(str(acos_b)+",")
        solar_angles_file.write(str(acos_c)+",")
        solar_angles_file.write(str(acos_d)+",")
        solar_angles_file.write(str(acos_e)+"\n")
        
        i=i+1
        
    solar_angles_file.close()

def solar_radiation(solar_angles_file, global_horizontal_radiation, direct_normal, diffuse_horizontal, DBT):
    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\solar_radiation_south.txt")
    radiation_soutn_file=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\solar_radiation_south.txt")
    next(solar_angles_file)
    
    for line in solar_angles_file:
        row_heading=line.split(",")
        incident_angle_south.append(row_heading[2])
        incident_angle_west.append(row_heading[3])
        incident_angle_east.append(row_heading[4])
        incident_angle_north.append(row_heading[5])
        incident_angle_horizontal.append(row_heading[6])
        
    incident_radiation_south_value=0
    incident_radiation_west_value=0
    incident_radiation_east_value=0
    incident_radiation_north_value=0
    incident_radiation_horizontal_value=0
    
    i=0
    for items in incident_angle_south:
        
        if math.cos(math.radians(float(incident_angle_south[i])))>0:

            incident_radiation_south_value = (math.cos(math.radians(float(incident_angle_south[i])))*float(direct_normal[i]))+(0.5*float(diffuse_horizontal[i]))+(0.5*0.2*float(global_horizontal_radiation[i]))

        else:
            
            incident_radiation_south_value =(0.5*float(diffuse_horizontal[i]))+(0.5*0.2*float(global_horizontal_radiation[i]))
            
        incident_radiation_south_value = round(incident_radiation_south_value,0)
            
        solar_radiation_south.append(incident_radiation_south_value)
        
        i=i+1
        
        
    i=0
    for items in incident_angle_west:
        
        if math.cos(math.radians(float(incident_angle_west[i])))>0:

            a=(0.5*float(diffuse_horizontal[i]))+(0.5*0.2*float(global_horizontal_radiation[i]))
            b=(math.cos(math.radians(float(incident_angle_west[i]))))*float(direct_normal[i])

            incident_radiation_west_value = (math.cos(math.radians(float(incident_angle_west[i])))*float(direct_normal[i]))+(0.5*float(diffuse_horizontal[i]))+(0.5*0.2*float(global_horizontal_radiation[i]))

        else:
            
            incident_radiation_west_value =(0.5*float(diffuse_horizontal[i]))+(0.5*0.2*float(global_horizontal_radiation[i]))
        
        incident_radiation_west_value = round (incident_radiation_west_value,0)   
        solar_radiation_west.append(incident_radiation_west_value)
        radiation_soutn_file.write(str(incident_radiation_west_value)+"\n")
        i=i+1
        
        
    i=0
    for items in incident_angle_east:
        
        if math.cos(math.radians(float(incident_angle_east[i])))>90:
            
            incident_radiation_east_value = math.cos(math.radians(float(incident_angle_east[i])))*float(direct_normal[i])+(0.5*float(diffuse_horizontal[i]))+(0.5*0.2*float(global_horizontal_radiation[i]))

        else:
            
            incident_radiation_east_value=(0.5*float(diffuse_horizontal[i]))+(0.5*0.2*float(global_horizontal_radiation[i]))
            
        solar_radiation_east.append(incident_radiation_east_value)
        incident_radiation_east_value = round(incident_radiation_east_value,0)
        
        i=i+1
             
    i=0
    for items in incident_angle_north:
        
        if math.cos(math.radians(float(incident_angle_north[i])))>90:
            
            incident_radiation_north_value = math.cos(math.radians(float(incident_angle_north[i])))*float(direct_normal[i])+(0.5*float(diffuse_horizontal[i]))+(0.5*0.2*float(global_horizontal_radiation[i]))

        else:
            
            incident_radiation_north_value =(0.5*float(diffuse_horizontal[i]))+(0.5*0.2*float(global_horizontal_radiation[i]))

            incident_radiation_north_value = round(incident_radiation_north_value,0)
            
        solar_radiation_north.append(incident_radiation_north_value)
        i=i+1
        
        
    i=0
    for items in incident_angle_horizontal:
        
        if math.cos(math.radians(float(incident_angle_horizontal[i])))>0:
            
            incident_radiation_horizontal_value = math.cos(math.radians(float(incident_angle_horizontal[i])))*float(direct_normal[i])+(float(diffuse_horizontal[i]))

        else:
            
            incident_radiation_horizontal_value =(float(diffuse_horizontal[i]))
            
        solar_radiation_horizontal.append(incident_radiation_horizontal_value)
        i=i+1
        
    
    return solar_radiation_south, solar_radiation_east, solar_radiation_west, solar_radiation_north, solar_radiation_horizontal

def envelope_gain(heat_carrying_capacity_of_air, building_area, building_height, set_point, Wall_area_south, Wall_area_west, Wall_area_east, Wall_area_north, Roof_area, wall_surf_refl, roof_surf_refl, wall_U_value, roof_U_value, DBT, solar_radiation_south, solar_radiation_east, solar_radiation_west, solar_radiation_north, solar_radiation_horizontal):

    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\envelope_gains_south.txt")
    envelope_gain_south_file=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\envelope_gains_south.txt")
    
    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\envelope_gains_west.txt")
    envelope_gain_west_file=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\envelope_gains_west.txt")
    
    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\envelope_gains_east.txt")
    envelope_gain_east_file=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\envelope_gains_east.txt")
    
    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\envelope_gains_north.txt")
    envelope_gain_north_file=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\envelope_gains_north.txt")
    
    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\envelope_gains_roof.txt")
    envelope_gain_roof_file=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\envelope_gains_roof.txt")
    
    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\envelope_gains.txt")
    envelope_gain_file=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\envelope_gains.txt")
    
    envelope_gain_file.write("temperature_rise_BAU,temperature_rise_ECBC\n")
    comfort_model_file=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\comfort_model.txt")
    
        
    i=0
    for items in solar_radiation_south:
        
        sol_air_temp_south_value=float(DBT[i])+(float(items)*float(wall_surf_refl))/22.7
        sol_air_temp_south.append(sol_air_temp_south_value)
        i=i+1
    
    i=0
    for items in solar_radiation_west:
        
        sol_air_temp_west_value=float(DBT[i])+(float(items)*float(wall_surf_refl))/22.7
        sol_air_temp_west.append(sol_air_temp_west_value)
        i=i+1
    
    i=0
    for items in solar_radiation_east:
        
        sol_air_temp_east_value=float(DBT[i])+(float(items)*float(wall_surf_refl))/22.7
        sol_air_temp_east.append(sol_air_temp_east_value)
        i=i+1
        
    i=0   
    for items in solar_radiation_north:
        
        sol_air_temp_north_value=float(DBT[i])+(float(items)*float(wall_surf_refl))/22.7
        sol_air_temp_north.append(sol_air_temp_north_value)
        i=i+1
        
    i=0   
    for items in solar_radiation_horizontal:
        
        sol_air_temp_roof_value=float(DBT[i])+((float(items)*float(roof_surf_refl))-(0.9*63))/22.7
        sol_air_temp_roof.append(sol_air_temp_roof_value)
        i=i+1


    i=0
    for items in sol_air_temp_south:
        
        if float(solar_radiation_south[i])>0:
            
            heat_gain_south_bau_value=(((float(Wall_area_south)*1.772)+((0.3*float(building_area)*float(building_height))/3))/1000)*((float(sol_air_temp_south[i])-float(set_point[i])))
            
            heat_gain_south_ecbc_value=(((float(Wall_area_south)*float(wall_U_value))+((0.3*float(building_area)*float(building_height))/3))/1000)*((float(sol_air_temp_south[i])-float(set_point[i])))
              
        else:
            
            heat_gain_south_bau_value = 0
            heat_gain_south_ecbc_value = 0

        heat_gain_south_bau.append(heat_gain_south_bau_value)
        heat_gain_south_ecbc.append(heat_gain_south_ecbc_value)
        
        heat_gain_south_bau_value= round(heat_gain_south_bau_value,2)
        heat_gain_south_ecbc_value = round(heat_gain_south_ecbc_value,2)
        
        envelope_gain_south_file.write(str(heat_gain_south_bau_value)+",")
        envelope_gain_south_file.write(str(heat_gain_south_ecbc_value)+"\n")
        i=i+1
        
    i=0
    for items in sol_air_temp_west:
        
        if float(solar_radiation_west[i])>0:
            
            heat_gain_west_bau_value=((float(Wall_area_west)*1.772)/1000)*(float(sol_air_temp_west[i])-float(set_point[i]))
            
            heat_gain_west_ecbc_value=(((float(Wall_area_west)*float(wall_U_value)))/1000)*(float(sol_air_temp_west[i])-float(set_point[i]))
            
        else:
            
            heat_gain_west_bau_value=0
            heat_gain_west_ecbc_value=0
            
        heat_gain_west_bau.append(heat_gain_west_bau_value)
        heat_gain_west_ecbc.append(heat_gain_west_ecbc_value)
        
        heat_gain_west_bau_value = round(heat_gain_west_bau_value,1)
        heat_gain_west_ecbc_value = round(heat_gain_west_ecbc_value,1)
        
        envelope_gain_west_file.write(str(heat_gain_west_bau_value)+",")
        envelope_gain_west_file.write(str(heat_gain_west_ecbc_value)+"\n")
        i=i+1
        
    i=0
    for items in sol_air_temp_east:

        if float(solar_radiation_east[i])>0:
        
            heat_gain_east_bau_value=(((float(Wall_area_east)*1.772))/1000)*(float(sol_air_temp_east[i])-float(set_point[i]))
            
            heat_gain_east_ecbc_value=(((float(Wall_area_east)*float(wall_U_value)))/1000)*(float(sol_air_temp_east[i])-float(set_point[i]))
            
        else:
            heat_gain_east_bau_value=0
            heat_gain_east_ecbc_value=0
            
        heat_gain_east_bau.append(heat_gain_east_bau_value)
        heat_gain_east_ecbc.append(heat_gain_east_ecbc_value)
        
        heat_gain_east_bau_value = round(heat_gain_east_bau_value,1)
        heat_gain_east_ecbc_value = round(heat_gain_east_ecbc_value,1)
        
        envelope_gain_east_file.write(str(heat_gain_east_bau_value)+",")
        envelope_gain_east_file.write(str(heat_gain_east_ecbc_value)+"\n")
        i=i+1
        
    i=0
    for items in sol_air_temp_north:

        if float(solar_radiation_north[i])>0:
        
            heat_gain_north_bau_value=(((float(Wall_area_north)*1.772))/1000)*(float(sol_air_temp_north[i])-float(set_point[i]))
            
            heat_gain_north_ecbc_value=((float(Wall_area_north)*float(wall_U_value))/1000)*(float(sol_air_temp_north[i])-float(set_point[i]))
            

        else:
            heat_gain_north_bau_value=0
            heat_gain_north_ecbc_value=0
            
        heat_gain_north_bau.append(heat_gain_north_bau_value)
        heat_gain_north_ecbc.append(heat_gain_north_ecbc_value)
        
        heat_gain_north_bau_value = round(heat_gain_north_bau_value,1)
        heat_gain_north_ecbc_value = round(heat_gain_north_ecbc_value,1)
        
        envelope_gain_north_file.write(str(heat_gain_north_bau_value)+",")
        envelope_gain_north_file.write(str(heat_gain_north_ecbc_value)+"\n")
        i=i+1
        
    i=0
    for items in sol_air_temp_roof:

        if float(solar_radiation_horizontal[i])>0:
        
            heat_gain_roof_bau_value=(((float(Roof_area)*2.942))/1000)*(float(sol_air_temp_roof[i])-float(set_point[i]))
            
            heat_gain_roof_ecbc_value=((float(Roof_area)*float(roof_U_value))/1000)*(float(sol_air_temp_roof[i])-float(set_point[i]))
            
        else:
            
            heat_gain_roof_bau_value=0
            heat_gain_roof_ecbc_value=0
            
        heat_gain_roof_bau.append(heat_gain_roof_bau_value)
        heat_gain_roof_ecbc.append(heat_gain_roof_ecbc_value)
        
        heat_gain_roof_bau_value = round(heat_gain_roof_bau_value,1)
        heat_gain_roof_ecbc_value = round(heat_gain_roof_ecbc_value,1)
        
        envelope_gain_roof_file.write(str(heat_gain_roof_bau_value)+",")
        envelope_gain_roof_file.write(str(heat_gain_roof_ecbc_value)+"\n")
        i=i+1
        
    i=0
    for items in heat_gain_south_bau:
        
        bau_envelope_value=heat_gain_south_bau[i]+heat_gain_west_bau[i]+heat_gain_east_bau[i]+heat_gain_north_bau[i]+heat_gain_roof_bau[i]
        envelope_gain_bau.append(bau_envelope_value)
        
        ecbc_envelope_value=heat_gain_south_ecbc[i]+heat_gain_west_ecbc[i]+heat_gain_east_ecbc[i]+heat_gain_north_ecbc[i]+heat_gain_roof_ecbc[i]
        envelope_gain_ecbc.append(ecbc_envelope_value)
        
        temp_rise_bau_value = bau_envelope_value/heat_carrying_capacity_of_air
        temp_rise_bau_value = round(temp_rise_bau_value, 2)
        envelope_gain_file.write(str(temp_rise_bau_value)+",")
        
        temp_rise_ecbc_value = ecbc_envelope_value/heat_carrying_capacity_of_air
        temp_rise_ecbc_value = round(temp_rise_ecbc_value, 2)
        envelope_gain_file.write(str(temp_rise_ecbc_value)+"\n")
        
        i=i+1
    envelope_gain_south_file.close()
    envelope_gain_west_file.close()
    envelope_gain_east_file.close()
    envelope_gain_north_file.close()
    envelope_gain_file.close()
    comfort_model_file.close()
    
def internal_load(heat_carrying_capacity_of_air, file_internal, occupancy, building_area, EPD, LPD, window_SHGC_nonnorth, window_SHGC_north, Win_area_south, Win_area_west, Win_area_east, Win_area_north, solar_radiation_south, solar_radiation_east, solar_radiation_west, solar_radiation_north, solar_radiation_horizontal):
    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\internal_gains.csv")
    internal_gains_file=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\internal_gains.csv")
    internal_gains_file.write("temperature_rise_internal_BAU,temperature_rise_internal_ECBC\n")

    next(file_internal)

    
    for line in file_internal:
        
        row_heading=line.split(",")
        lighting_schedule.append(row_heading[5])
        equipment_schedule.append(row_heading[4])

    equipment_energy=0
    lighting_energy=0
    i=0
    for items in occupancy_schedule:
        
        occu_heat_gain_value=float(occupancy_schedule[i])*70*int(occupancy)/1000
        occu_heat_gain.append(occu_heat_gain_value)
        
        equip_heat_gain_value=float(equipment_schedule[i])*float(EPD)*int(building_area)/1000
        equip_heat_gain.append(equip_heat_gain_value)
        equipment_energy=equipment_energy+equip_heat_gain_value
        
        light_heat_gain_value=float(lighting_schedule[i])*float(LPD)*int(building_area)/1000
        light_heat_gain.append(light_heat_gain_value)
        lighting_energy=lighting_energy+light_heat_gain_value
        i=i+1
    equipment_energy = equipment_energy / float(building_area)
    lighting_energy = lighting_energy / float(building_area)

    i=0    
    for items in solar_radiation_south:
        
        radiant_gain_south_BAU_value=float(solar_radiation_south[i])*float(Win_area_south)*0.82/1000
        radiant_gain_south_BAU.append(radiant_gain_south_BAU_value)
        
        radiant_gain_south_ECBC_value=float(solar_radiation_south[i])*float(Win_area_south)*float(window_SHGC_nonnorth)/1000
        radiant_gain_south_ECBC.append(radiant_gain_south_ECBC_value)
        
        i=i+1
    
    i=0    
    for items in solar_radiation_west:
        
        radiant_gain_west_BAU_value=float(solar_radiation_west[i])*float(Win_area_west)*0.82/1000
        radiant_gain_west_BAU.append(radiant_gain_west_BAU_value)
        
        radiant_gain_west_ECBC_value=float(solar_radiation_west[i])*float(Win_area_west)*float(window_SHGC_nonnorth)/1000
        radiant_gain_west_ECBC.append(radiant_gain_west_ECBC_value)
        
        i=i+1

    i=0    
    for items in solar_radiation_east:
        
        radiant_gain_east_BAU_value=float(solar_radiation_east[i])*float(Win_area_east)*0.82/1000
        radiant_gain_east_BAU.append(radiant_gain_east_BAU_value)
        
        radiant_gain_east_ECBC_value=float(solar_radiation_east[i])*float(Win_area_east)*float(window_SHGC_nonnorth)/1000
        radiant_gain_east_ECBC.append(radiant_gain_east_ECBC_value)
        
        i=i+1

    i=0    
    for items in solar_radiation_north:
        
        radiant_gain_north_BAU_value=float(solar_radiation_north[i])*float(Win_area_north)*0.82/1000
        radiant_gain_north_BAU.append(radiant_gain_north_BAU_value)
        
        radiant_gain_north_ECBC_value=float(solar_radiation_north[i])*float(Win_area_north)*float(window_SHGC_north)/1000
        radiant_gain_north_ECBC.append(radiant_gain_north_ECBC_value)
        i=i+1
        
    internal_temp_rise_bau_value=0
    internal_temp_rise_ecbc_value=0
    i=0
    for items in occu_heat_gain:
        
        total_internal_gain_bau_value=(occu_heat_gain[i]+equip_heat_gain[i]+light_heat_gain[i])+radiant_gain_south_BAU[i]+radiant_gain_west_BAU[i]+radiant_gain_east_BAU[i]+radiant_gain_north_BAU[i]
        total_internal_gain_bau_value = round(total_internal_gain_bau_value,2)
        total_internal_gain_bau.append(total_internal_gain_bau_value)
        internal_gains_file.write(str(total_internal_gain_bau_value)+",")
        
        total_internal_gain_ecbc_value=occu_heat_gain[i]+equip_heat_gain[i]+light_heat_gain[i]+radiant_gain_south_ECBC[i]+radiant_gain_west_ECBC[i]+radiant_gain_east_ECBC[i]+radiant_gain_north_ECBC[i]
        total_internal_gain_ecbc_value = round(total_internal_gain_ecbc_value,2)
        total_internal_gain_ecbc.append(total_internal_gain_ecbc_value)
        internal_gains_file.write(str(total_internal_gain_ecbc_value)+",")
        
        internal_temp_rise_bau_value=float(total_internal_gain_bau_value)/float(heat_carrying_capacity_of_air)
        internal_temp_rise_bau_value = round(internal_temp_rise_bau_value,2)
        internal_temp_rise_bau.append(internal_temp_rise_bau_value)
        internal_gains_file.write(str(internal_temp_rise_bau_value)+",")
        
        internal_temp_rise_ecbc_value=float(total_internal_gain_ecbc_value)/float(heat_carrying_capacity_of_air)
        internal_temp_rise_ecbc_value = round(internal_temp_rise_ecbc_value,2)
        internal_temp_rise_ecbc.append(internal_temp_rise_ecbc_value)
        internal_gains_file.write(str(internal_temp_rise_ecbc_value)+"\n")
        i=i+1
        
    internal_gains_file.close()
    return occu_heat_gain, equip_heat_gain, light_heat_gain,  total_internal_gain_bau, total_internal_gain_ecbc, internal_temp_rise_bau_value, internal_temp_rise_ecbc, equipment_energy, lighting_energy
    
    
def mitigation_at_night(file_internal, time_constant_BAU, time_constant_ECBC, thermal_capacitance, heat_carrying_capacity_of_air, occupancy_schedule, set_point):
    
    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\mitigation_at_night.txt")
    mitigation_at_night_file=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\mitigation_at_night.txt")
    
    comfort_model_file=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\comfort_model.txt")
    
    avg_night_daily_file=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\avg_night_daily.txt")
    
    file_internal=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\Internal_load.csv")
    
    next(file_internal)
        
    for items in avg_night_daily_file:
        avg_night_daily.append(items)
        
    i=0
    for items in occupancy_schedule:

        a=-math.exp((-1*12/float(time_constant_BAU))-1)
        b=(float(set_point[i])-float(avg_night_daily[i]))
        mitigation_BAU_value=((a*b)*float(thermal_capacitance))/(24*3600*float(heat_carrying_capacity_of_air))
        mitigation_BAU_value = round(mitigation_BAU_value,2)
          
        mitigation_ECBC_value = -math.exp((-1*12/float(time_constant_ECBC))-1)*(float(set_point[i])-float(avg_night_daily[i]))*float(thermal_capacitance)/(24*3600*float(heat_carrying_capacity_of_air))
        mitigation_ECBC_value = round(mitigation_ECBC_value,2)
        
        mitigation_at_night_file.write(str(mitigation_BAU_value)+",")
        mitigation_at_night_BAU.append(mitigation_BAU_value)
        
        mitigation_at_night_file.write(str(mitigation_ECBC_value)+"\n")
        mitigation_at_night_ECBC.append(mitigation_ECBC_value)
        i=i+1
        
    return mitigation_at_night_BAU, mitigation_at_night_ECBC

    mitigation_at_night_file.close()
    file_internal.close()
    avg_night_daily_file.close()
    comfort_model_file.close()

def fan_gain(heat_carrying_capacity_of_air):

    fan_temp_rise=(float(heat_carrying_capacity_of_air)/1.224)*1.5/(float(heat_carrying_capacity_of_air)*0.7)
    fan_temp_rise = round(fan_temp_rise, 2)

    return fan_temp_rise

def latent_load(occupancy, heat_carrying_capacity_of_air, Out_humdrat, In_humdrat, fan_temp_rise, atm_prs, set_point, occupancy_schedule, total_internal_gain_bau, total_internal_gain_ecbc, internal_temp_rise_bau_value, internal_temp_rise_ecbc, mitigation_at_night_BAU, mitigation_at_night_ECBC):
    
    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\latent_gains_bau.txt")
    latent_temp_rise_bau_file=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\latent_gains_bau.txt")

    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\latent_gains_ecbc.txt")
    latent_temp_rise_ecbc_file=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\latent_gains_ecbc.txt")
    
    #file_internal=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\Internal_load.csv")
    #next(file_internal)

    internal_gains_file=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\internal_gains.csv")
    next(internal_gains_file)

    envelope_gain_file=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\envelope_gains.txt")
    next(envelope_gain_file)

    mitigation_at_night_file=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\mitigation_at_night.txt")
    

    for line in envelope_gain_file:
        row_heading=line.split(",")
        envelope_temperature_rise_BAU.append(row_heading[0])
        envelope_temperature_rise_ECBC.append(row_heading[1])

    #for line in mitigation_at_night_file:
        #row_heading=line.split(",")
        #mitigation_at_night_BAU.append(row_heading[0])
        #mitigation_at_night_ECBC.append(row_heading[1])
        
    i=0
    for items in occupancy_schedule:
        latent_load_occupant_value=46.5*float(occupancy_schedule[i])*float(occupancy)/1000
        occupant_moisture_gain_value=latent_load_occupant_value/(float(heat_carrying_capacity_of_air)*2450)
        occupant_moisture_gain.append(occupant_moisture_gain_value)
        
        i=i+1
        
    i=0
    for items in set_point:

        supply_air_temp_BAU_value=float(set_point[i])-float(envelope_temperature_rise_BAU[i])-float(internal_temp_rise_bau[i])-float(mitigation_at_night_BAU[i])-float(fan_temp_rise) 
        supply_air_temp_BAU.append(supply_air_temp_BAU_value)
        
        supply_air_temp_ECBC_value=float(set_point[i])-float(envelope_temperature_rise_ECBC[i])-float(internal_temp_rise_ecbc[i])-float(mitigation_at_night_ECBC[i])-float(fan_temp_rise)
        supply_air_temp_ECBC.append(supply_air_temp_ECBC_value)

        C1 = -5674.5359
    
        C2 = 6.3925247
    
        C3 = -0.009677843
    
        C4 = 0.00000062215701
    
        C5 = 2.0747825E-09
    
        C6 = -9.484024E-13
    
        C7 = 4.1635019
    
        C8 = -5800.2206
    
        C9 = 1.3914993
    
        C10 = -0.048640239
    
        C11 = 0.000041764768
    
        C12 = -0.000000014452093
    
        C13 = 6.5459673
        
        Tk_bau = supply_air_temp_BAU_value + 273.5
        Tk_ecbc=supply_air_temp_ECBC_value + 273.5
        
        if Tk_bau <= 273.5:
            
            Pws_bau_value = math.exp(C1/Tk_bau + C2 + C3*Tk_bau + C4*Tk_bau**2 + C5*Tk_bau**3 + C6*Tk_bau**4 + C7*math.log(Tk_bau)) / 1000

        else:

            Pws_bau_value = math.exp(C8/Tk_bau + C9 + C10*Tk_bau + C11*Tk_bau**2 + C12*Tk_bau**3 + C13*math.log(Tk_bau)) / 1000
        
        if Tk_ecbc <= 273.5:
            
            Pws_ecbc_value = math.exp(C1/Tk_ecbc + C2 + C3*Tk_ecbc + C4*Tk_ecbc**2 + C5*Tk_ecbc**3 + C6*Tk_ecbc**4 + C7*math.log(Tk_ecbc)) / 1000

        else:

            Pws_ecbc_value = math.exp(C8/Tk_ecbc + C9 + C10*Tk_ecbc + C11*Tk_ecbc**2 + C12*Tk_ecbc**3 + C13*math.log(Tk_ecbc)) / 1000   
            
        supply_humdrat_bau_value=0.621945*0.6*float(Pws_bau_value)/((float(atm_prs[i])/1000)-0.6*float(Pws_bau_value))
        supply_humdrat_ecbc_value=0.621945*0.6*float(Pws_ecbc_value)/((float(atm_prs[i])/1000)-0.6*float(Pws_ecbc_value))
        
        
        supply_humdrat_bau.append(supply_humdrat_bau_value)
        supply_humdrat_ecbc.append(supply_humdrat_ecbc_value)
        
        i=i+1

    i=0    
    for items in Out_humdrat:
        
        infiltration_value = float(Out_humdrat[i])*0.1
        infiltration.append(infiltration_value)
        
        on_coil_moisture_bau_value = infiltration_value+occupant_moisture_gain[i]+supply_humdrat_bau[i]
        on_coil_moisture_bau.append(on_coil_moisture_bau_value)

        on_coil_moisture_ecbc_value = infiltration_value+occupant_moisture_gain[i]+supply_humdrat_ecbc[i]
        on_coil_moisture_ecbc.append(on_coil_moisture_ecbc_value)
       
        i=i+1
        
    C14 = 6.54

    C15 = 14.526

    C16 = 0.7389

    C17 = 0.09486

    C18 = 0.4569

    i=0
    for items in atm_prs:

        Pw_bau_value = (float(atm_prs[i])/1000)*float(on_coil_moisture_bau[i])/(0.621945+float(on_coil_moisture_bau[i]))
        Pw_ecbc_value = (float(atm_prs[i])/1000)*float(on_coil_moisture_ecbc[i])/(0.621945+float(on_coil_moisture_ecbc[i]))

        alpha_bau_value = math.log(Pw_bau_value)
        alpha_ecbc_value = math.log(Pw_ecbc_value)

        Tdp1_bau_value = C14 + C15*alpha_bau_value + C16*alpha_bau_value**2 + C17*alpha_bau_value**3 + C18*Pw_bau_value**0.1984

        Tdp2_bau_value = 6.09 + 12.608*alpha_bau_value + 0.4959*alpha_bau_value**2

        Tdp1_ecbc_value = C14 + C15*alpha_ecbc_value + C16*alpha_ecbc_value**2 + C17*alpha_ecbc_value**3 + C18*Pw_ecbc_value**0.1984

        Tdp2_ecbc_value = 6.09 + 12.608*alpha_ecbc_value + 0.4959*alpha_ecbc_value**2

        if Tdp1_bau_value >= 0:

            dew_point_bau_value = Tdp1_bau_value

        else:
            
            dew_point_bau_value = Tdp2_bau_value

            

        if Tdp1_ecbc_value >= 0:

            dew_point_ecbc_value = Tdp1_ecbc_value

        else:

            dew_point_ecbc_value = Tdp2_ecbc_value

            
        dew_point_bau_value = round(dew_point_bau_value,2)
        dew_point_ecbc_value = round(dew_point_ecbc_value,2)
        Room_air_dew_point_bau.append(dew_point_bau_value)
        Room_air_dew_point_ecbc.append(dew_point_ecbc_value)
            
        i=i+1

    i=0
    for items in Room_air_dew_point_bau:
        if float(items)<float(set_point[i]):
            
            latent_temp_rise_bau_value = 2400*(on_coil_moisture_bau[i] - supply_humdrat_bau[i])

        else:

            latent_temp_rise_bau_value = 0
        latent_temp_rise_bau_value = round(latent_temp_rise_bau_value,2)
        latent_temp_rise_bau.append(latent_temp_rise_bau_value)
        latent_temp_rise_bau_file.write(str(latent_temp_rise_bau_value)+"\n")
        i=i+1
        
    i=0
    for items in Room_air_dew_point_ecbc:
        if float(items)<float(set_point[i]):
            
            latent_temp_rise_ecbc_value = 2400*(on_coil_moisture_ecbc[i] - supply_humdrat_ecbc[i])

        else:

            latent_temp_rise_ecbc_value = 0
        latent_temp_rise_ecbc_value = round(latent_temp_rise_ecbc_value,2)
        latent_temp_rise_ecbc.append(latent_temp_rise_ecbc_value)
        latent_temp_rise_ecbc_file.write(str(latent_temp_rise_ecbc_value)+"\n")
        i=i+1

    return latent_temp_rise_bau, latent_temp_rise_ecbc
    latent_temp_rise_bau_file.close()
    latent_temp_rise_ecbc_file.close()
    #file_internal.close()
    internal_gains_file.close()
    envelope_gain_file.close()
    mitigation_at_night_file.close()

def CDD_std(DBT):
    
    CDD_18_3=0
    
    i=0
    for items in DBT:
        
        CDD_value=float(DBT[i])-18.3
        
        if CDD_value>0:
            CDD_18_3 = (CDD_18_3 + CDD_value)
        i=i+1
    
    CDD_18_3 = CDD_18_3 / 24
    CDD_18_3 = round(CDD_18_3,0)
    return CDD_18_3
    

def base_temperature(heat_carrying_capacity_of_air, occupancy_schedule, DBT, set_point, envelope_temperature_rise_BAU, envelope_temperature_rise_ECBC, latent_temp_rise_bau, latent_temp_rise_ecbc, mitigation_at_night_BAU, mitigation_at_night_ECBC, internal_temp_rise_bau, internal_temp_rise_ecbc, fan_temp_rise, building_area):

    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\BPT.csv")
    bpt_file= open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\BPT.csv")
    bpt_file.write("bpt_bau,bpt_ecbc\n")
    CDD_AC_BAU = 0
    CDD_AC_ECBC = 0
    i=0
    for items in set_point:
        bpt_bau_value = float(set_point[i]) - float(envelope_temperature_rise_BAU[i]) - float(latent_temp_rise_bau[i]) - float(mitigation_at_night_BAU[i])- float(internal_temp_rise_bau[i]) - float(fan_temp_rise)
        bpt_ecbc_value = float(set_point[i]) - float(envelope_temperature_rise_ECBC[i]) - float(latent_temp_rise_ecbc[i]) - float(mitigation_at_night_ECBC[i])- float(internal_temp_rise_ecbc[i]) - float(fan_temp_rise)

        bpt_file.write(str(bpt_bau_value)+",")
        bpt_file.write(str(bpt_ecbc_value)+"\n")
        bpt_bau.append(bpt_bau_value)
        bpt_ecbc.append(bpt_ecbc_value)

        CDH_bau_value = float(DBT[i]) - bpt_bau_value
        CDH_bau.append(CDH_bau_value)
        CDH_ECBC_value = float(DBT[i]) - bpt_ecbc_value
        CDH_ECBC.append(CDH_ECBC_value)
        #print("CDH_BAU: "+str(CDH_bau))
        #print("CDH_ECBC: "+str(CDH_ECBC))
        
        if CDH_bau_value>0:
            
            CDD_AC_BAU = CDD_AC_BAU + CDH_bau_value

        if CDH_ECBC_value>0:
            
            CDD_AC_ECBC = CDD_AC_ECBC + CDH_ECBC_value
    
        i=i+1
    CDD_AC_BAU=CDD_AC_BAU/24
    CDD_AC_BAU = round(CDD_AC_BAU,0)
    
    CDD_AC_ECBC=CDD_AC_ECBC/24
    CDD_AC_ECBC = round(CDD_AC_ECBC,0)
    
    CDD_BAU_energy=0
    CDD_ECBC_energy=0
    i=0
    for items in occupancy_schedule:
        
        if float(occupancy_schedule[i])>0:
            
            CDD_BAU_energy = CDD_BAU_energy + CDH_bau[i]
            CDD_ECBC_energy = CDD_ECBC_energy + CDH_ECBC[i]
        i=i+1

    cooling_energy_AC_bau= float(CDD_BAU_energy)*float(heat_carrying_capacity_of_air)/3.1
    cooling_energy_AC_bau = round(cooling_energy_AC_bau,0)
    cooling_energy_AC_bau = cooling_energy_AC_bau / float(building_area)
    
    cooling_energy_AC_ecbc= float(CDD_ECBC_energy)*float(heat_carrying_capacity_of_air)/3.1
    cooling_energy_AC_ecbc = round(cooling_energy_AC_ecbc,0)
    cooling_energy_AC_ecbc = cooling_energy_AC_ecbc / float(building_area)

    #print(cooling_energy_AC_bau)
    #print(cooling_energy_AC_ecbc)
    bpt_file.close()
    return bpt_bau, bpt_ecbc, CDD_BAU_energy, CDD_ECBC_energy, cooling_energy_AC_bau, cooling_energy_AC_ecbc, CDD_AC_BAU, CDD_AC_ECBC

def comfort_cooling(occupancy_schedule, DBT, wind_speed, bpt_bau, bpt_ecbc, heat_carrying_capacity_of_air, building_area):


    RCDD_comfcool_bau = 0
    RCDD_comfcool_ecbc = 0
    i=0
    for items in DBT:

        indoor_wind_speed_value=0.3*float(wind_speed[i])
        temp_dec_perc_value=(2.319*indoor_wind_speed_value)+0.4816
        indoor_temp_perc_value=float(DBT[i])-temp_dec_perc_value
        
        RCDH_comfcool_bau_value = indoor_temp_perc_value - bpt_bau[i]
        RCDH_comfcool_ecbc_value = indoor_temp_perc_value - bpt_ecbc[i]

        RCDH_comfcool_bau.append(RCDH_comfcool_bau_value)
        RCDH_comfcool_ecbc.append(RCDH_comfcool_ecbc_value)

        if RCDH_comfcool_bau_value > 0:
            RCDD_comfcool_bau = RCDD_comfcool_bau + RCDH_comfcool_bau_value

        if RCDH_comfcool_ecbc_value > 0:
            RCDD_comfcool_ecbc = RCDD_comfcool_ecbc + RCDH_comfcool_ecbc_value
        i=i+1
        
    RCDD_comfcool_bau = RCDD_comfcool_bau/24
    RCDD_comfcool_bau = round(RCDD_comfcool_bau,0)
    
    RCDD_comfcool_ecbc = RCDD_comfcool_ecbc/24
    RCDD_comfcool_ecbc = round(RCDD_comfcool_ecbc,0)
        
    RCDD_comfcool_BAU_energy=0
    RCDD_comfcool_ECBC_energy=0
    i=0
    for items in occupancy_schedule:
        
        if float(occupancy_schedule[i])>0:
            
            RCDD_comfcool_BAU_energy = RCDD_comfcool_BAU_energy + RCDH_comfcool_bau[i]
            RCDD_comfcool_ECBC_energy = RCDD_comfcool_ECBC_energy + RCDH_comfcool_ecbc[i]
        i=i+1
        
    cooling_energy_comfcool_bau=float(RCDD_comfcool_BAU_energy)*float(heat_carrying_capacity_of_air)/3.1
    cooling_energy_comfcool_bau = round(cooling_energy_comfcool_bau,0)
    cooling_energy_comfcool_bau = cooling_energy_comfcool_bau / float(building_area)
    
    cooling_energy_comfcool_ecbc=float(RCDD_comfcool_ECBC_energy)*float(heat_carrying_capacity_of_air)/3.1
    cooling_energy_comfcool_ecbc = round(cooling_energy_comfcool_ecbc,0)
    cooling_energy_comfcool_ecbc = cooling_energy_comfcool_ecbc / float(building_area)
    return RCDH_comfcool_bau, RCDH_comfcool_ecbc, cooling_energy_comfcool_bau, cooling_energy_comfcool_ecbc, RCDD_comfcool_bau, RCDD_comfcool_ecbc


def evaporative_cooling(DBT, set_point, bpt_bau, bpt_ecbc, Tw, occupancy_schedule, heat_carrying_capacity_of_air, building_area, evap_power):
    i=0
    for items in occupancy_schedule:

        WBD_value = float(DBT[i]) - float(Tw[i])
        
        if float(occupancy_schedule[i]) > 0:
            evap_temp_value = float(DBT[i])-(WBD_value*0.9)
        else:
            evap_temp_value = float(DBT[i])
        evap_temp.append(evap_temp_value)

        
            
        if evap_temp[i] > float(bpt_bau[i]):
            RCDH_precool_bau_value = evap_temp[i] - float(bpt_bau[i])
        else:
            RCDH_precool_bau_value = 0

        RCDH_precool_bau.append(RCDH_precool_bau_value)

        if evap_temp[i] > float(bpt_ecbc[i]):
            RCDH_precool_ecbc_value = evap_temp[i] - float(bpt_ecbc[i])
        else:
            RCDH_precool_ecbc_value = 0

        RCDH_precool_ecbc.append(RCDH_precool_ecbc_value)

        
        
        if evap_temp[i] > float(bpt_bau[i]):
            RCDH_switch_bau_value = float(DBT[i]) - float(bpt_bau[i])
        else:
            RCDH_switch_bau_value = 0

        RCDH_switch_bau.append(RCDH_switch_bau_value)

        if evap_temp[i] > float(bpt_ecbc[i]):
            RCDH_switch_ecbc_value = float(DBT[i]) - float(bpt_ecbc[i])
        else:
            RCDH_switch_ecbc_value = 0

        RCDH_switch_ecbc.append(RCDH_switch_ecbc_value)
        
        i=i+1

    pre_cool_RCDD_bau=0
    i=0
    for items in RCDH_precool_bau:
        if float(RCDH_precool_bau[i])>0:
            pre_cool_RCDD_bau = pre_cool_RCDD_bau + RCDH_precool_bau[i]
        i=i+1
    pre_cool_RCDD_bau = pre_cool_RCDD_bau / 24
    
    pre_cool_RCDD_ecbc=0
    i=0
    for items in RCDH_precool_ecbc:
        if float(RCDH_precool_ecbc[i])>0:
            pre_cool_RCDD_ecbc=pre_cool_RCDD_ecbc + RCDH_precool_ecbc[i]
        i=i+1
    pre_cool_RCDD_ecbc = pre_cool_RCDD_ecbc / 24

    
    evap_cool_RCDD_ecbc=0
    i=0
    for items in RCDH_switch_ecbc:
        if float(RCDH_switch_ecbc[i])>0:
            evap_cool_RCDD_ecbc=evap_cool_RCDD_ecbc + RCDH_switch_ecbc[i]
        i=i+1
    evap_cool_RCDD_ecbc = evap_cool_RCDD_ecbc / 24

    
    evap_cool_RCDD_bau=0
    i=0
    for items in RCDH_switch_bau:
        if float(RCDH_switch_bau[i])>0:
            evap_cool_RCDD_bau = evap_cool_RCDD_bau + RCDH_switch_bau[i]
        i=i+1
    evap_cool_RCDD_bau = evap_cool_RCDD_bau / 24

    i=0
    for items in occupancy_schedule:
        if float(occupancy_schedule[i]) > 0:
            evap_precool_energy_bau_hourly_value =  RCDH_precool_bau[i] * heat_carrying_capacity_of_air / 3.1
        else:
            evap_precool_energy_bau_hourly_value = 0
        evap_precool_energy_hourly_bau.append(evap_precool_energy_bau_hourly_value)
        i=i+1
        
    evap_precool_energy_bau = 0
    i=0
    for items in evap_precool_energy_hourly_bau:
        if float(evap_precool_energy_hourly_bau[i]) > 0:
            evap_precool_energy_bau = evap_precool_energy_bau + evap_precool_energy_hourly_bau[i]
        i=i+1
    evap_precool_energy_bau = evap_precool_energy_bau / float(building_area)

    
    i=0
    for items in occupancy_schedule:
        if float(occupancy_schedule[i]) > 0:
            evap_precool_energy_hourly_ecbc_value =  RCDH_precool_ecbc[i] * heat_carrying_capacity_of_air / 3.1
        else:
            evap_precool_energy_hourly_ecbc_value = 0
        evap_precool_energy_hourly_ecbc.append(evap_precool_energy_hourly_ecbc_value)
        i=i+1

    evap_precool_energy_ecbc = 0
    i=0
    for items in evap_precool_energy_hourly_ecbc:
        
        if float(evap_precool_energy_hourly_ecbc[i]) > 0:
            evap_precool_energy_ecbc = evap_precool_energy_ecbc + evap_precool_energy_hourly_ecbc[i]
        i=i+1
    evap_precool_energy_ecbc = evap_precool_energy_ecbc / float(building_area)

        
    i=0
    for items in occupancy_schedule:
        if float(occupancy_schedule[i]) > 0:
            evap_cool_energy_hourly_bau_value = RCDH_switch_bau[i] * heat_carrying_capacity_of_air / 3.1
        else:
            evap_cool_energy_hourly_bau_value = 0
        evap_cool_energy_hourly_bau.append(evap_cool_energy_hourly_bau_value)
        i=i+1

    evap_cool_energy_bau = 0
    i=0
    for items in evap_cool_energy_hourly_bau:
        if float(evap_cool_energy_hourly_bau[i]) > 0:
            evap_cool_energy_bau = evap_cool_energy_bau + evap_cool_energy_hourly_bau[i]
        i=i+1
    evap_cool_energy_bau = evap_cool_energy_bau / float(building_area)

    i=0
    for items in occupancy_schedule:
        if float(occupancy_schedule[i]) > 0:
            evap_cool_energy_hourly_ecbc_value = RCDH_switch_ecbc[i] * heat_carrying_capacity_of_air / 3.1
        else:
            evap_cool_energy_hourly_ecbc_value = 0
        evap_cool_energy_hourly_ecbc.append(evap_cool_energy_hourly_ecbc_value)
        i=i+1

    evap_cool_energy_ecbc = 0
    i=0
    for items in evap_cool_energy_hourly_ecbc:
        if float(evap_cool_energy_hourly_ecbc[i]) > 0:
            evap_cool_energy_ecbc = evap_cool_energy_ecbc + evap_cool_energy_hourly_ecbc[i]
        i=i+1
    evap_cool_energy_ecbc = evap_cool_energy_ecbc / float(building_area)

    evap_precooling_energy = 0
    i=0
    for items in occupancy_schedule:
        if float(occupancy_schedule[i])>0:
            evap_precooling_energy_value = float(evap_power)
        else:
            evap_precooling_energy_value = 0

        evap_precooling_energy = evap_precooling_energy + evap_precooling_energy_value
        i=i+1
    evap_precooling_energy = evap_precooling_energy / float(building_area)

    print("evap_precooling_energy_: "+str(evap_precooling_energy))

    
    evap_cooling_energy_bau = 0
    i=0
    for items in  bpt_bau:
        if float(occupancy_schedule[i])>0:
            if float(evap_temp[i])<= float(bpt_bau[i]):
                evap_cooling_energy_bau_value = float(evap_power)
            else:
                evap_cooling_energy_bau_value = 0
        else:
            evap_cooling_energy_bau_value = 0

        evap_cooling_energy_bau = evap_cooling_energy_bau + evap_cooling_energy_bau_value
        i=i+1
    evap_cooling_energy_bau = evap_cooling_energy_bau / float(building_area)
    print("evap_cooling_energy_bau: "+str(evap_cooling_energy_bau))

    evap_cooling_energy_ecbc = 0
    i=0
    for items in  bpt_ecbc:
        if float(occupancy_schedule[i])>0:
            if float(evap_temp[i])<= float(bpt_ecbc[i]):
                evap_cooling_energy_ecbc_value = float(evap_power)
            else:
                evap_cooling_energy_ecbc_value = 0
        else:
            evap_cooling_energy_ecbc_value = 0

        evap_cooling_energy_ecbc = evap_cooling_energy_ecbc + evap_cooling_energy_ecbc_value
        i=i+1
    evap_cooling_energy_ecbc = evap_cooling_energy_ecbc / float(building_area)
    print("evap_cooling_energy_ecbc: "+str(evap_cooling_energy_ecbc))
    
          
    return pre_cool_RCDD_bau, pre_cool_RCDD_ecbc, evap_cool_RCDD_ecbc, evap_cool_RCDD_bau, evap_precool_energy_bau, evap_precool_energy_ecbc, evap_cool_energy_bau, evap_cool_energy_ecbc, evap_precooling_energy, evap_cooling_energy_bau, evap_cooling_energy_ecbc

def night_ventilation(DBT, abc, seas_avg, monthly_min, minimum_seasonal, daily_max, daily_min, avg_dbt, bpt_bau, bpt_ecbc, heat_carrying_capacity_of_air, building_area):
    
    i=0
    while i<len(avg_dbt):
        j=0
        while j<24:
            average_DBT.append(avg_dbt[i])
            j=j+1
        i=i+1

    i=0
    while i<len(daily_max):
        j=0
        while j<24:
            daily_maximum.append(daily_max[i])
            j=j+1
        i=i+1

    i=0
    while i<len(daily_min):
        j=0
        while j<24:
            daily_minimum.append(daily_min[i])
            j=j+1
        i=i+1
   
    i=0
    for items in seas_avg:
        TinmaxECBC_value = float(seas_avg[i])+6-(0.006*float(seas_avg[i]))+ (0.9*(float(average_DBT[i])-float(seas_avg[i])))
        TinminECBC_value = float(minimum_seasonal[i]) + 5.8 - (0.182 * float(minimum_seasonal[i])) + (0.75 * (float(daily_minimum[i]) - float(minimum_seasonal[i]))) + (0.26 * (float(average_DBT[i-24] - float(daily_minimum[i]))))
        TinmaxBAU_value = float(seas_avg[i]) + 3.9 + (0.108 * float(seas_avg[i])) + (0.8 * (float(average_DBT[i])-float(seas_avg[i])))
        TinminBAU_value = float(minimum_seasonal[i]) + 4.9 - (0.068 * float(minimum_seasonal[i])) + (0.75 * (float(daily_minimum[i]) - float(minimum_seasonal[i]))) + (0.23 * (float(average_DBT[i-24] - float(daily_minimum[i]))))
        amplitude_bau_value = TinmaxBAU_value - TinminBAU_value
        amplitude_ecbc_value = TinmaxECBC_value - TinminECBC_value
        TinmaxECBC.append(TinmaxECBC_value)
        TinminECBC.append(TinminECBC_value)
        TinminBAU.append(TinminBAU_value)
        TinmaxBAU.append(TinmaxBAU_value)
        amplitude_bau.append(amplitude_bau_value)
        amplitude_ecbc.append(amplitude_ecbc_value)
        i=i+1

    i=0
    while i<365:
        j=0
        while j<24:
            angle_degree_value = j*15
            angle_radian_value = math.radians(angle_degree_value)
            sine_theeta_value = math.sin(angle_radian_value)
            angle_radian.append(angle_radian_value)
            sine_theeta.append(sine_theeta_value)
            j=j+1
        i=i+1

    i=0
    for items in TinmaxBAU:
        Tin_hourly_bau_value = (amplitude_bau[i]/2*math.sin(-1*angle_radian[i]))+(TinmaxBAU[i]+TinminBAU[i])/2
        Tin_hourly_bau.append(Tin_hourly_bau_value)
        i=i+1

    i=0
    for items in TinmaxECBC:
        Tin_hourly_ecbc_value = (amplitude_ecbc[i]/2*math.sin(-1*angle_radian[i]))+(TinmaxECBC[i]+TinminECBC[i])/2
        Tin_hourly_ecbc.append(Tin_hourly_ecbc_value)
        i=i+1

    RCDD_night_bau = 0
    i=0
    for items in Tin_hourly_bau:
        if float(Tin_hourly_bau[i]) > float(set_point[i]):
            RCDH_night_bau_value = Tin_hourly_bau[i] - bpt_bau[i]
        else:
            RCDH_night_bau_value = 0
        RCDH_night_bau.append(RCDH_night_bau_value)

        if RCDH_night_bau_value > 0:
            RCDD_night_bau = RCDD_night_bau + RCDH_night_bau_value
        i=i+1

    nightcool_energy_bau = 0
    i=0
    for items in RCDH_night_bau:
        if float(occupancy_schedule[i]) > 0:
            nightcool_energy_hourly_bau_value = RCDH_night_bau[i] * heat_carrying_capacity_of_air / 3.1
        else:
            nightcool_energy_hourly_bau_value = 0
        
        nightcool_energy_hourly_bau.append(nightcool_energy_hourly_bau_value)
        nightcool_energy_bau = nightcool_energy_bau + nightcool_energy_hourly_bau_value
        i=i+1
    nightcool_energy_bau = nightcool_energy_bau / float(building_area)

    RCDD_night_ecbc = 0
    i=0
    for items in Tin_hourly_ecbc:
        if float(Tin_hourly_ecbc[i]) > float(set_point[i]):
            RCDH_night_ecbc_value = Tin_hourly_ecbc[i] - bpt_ecbc[i]
        else:
            RCDH_night_ecbc_value = 0
        RCDH_night_ecbc.append(RCDH_night_ecbc_value)

        if RCDH_night_ecbc_value > 0:
            RCDD_night_ecbc = RCDD_night_ecbc + RCDH_night_ecbc_value
        i=i+1

    nightcool_energy_ecbc = 0
    i=0
    for items in RCDH_night_ecbc:
        if float(occupancy_schedule[i]) > 0:
            nightcool_energy_hourly_ecbc_value = RCDH_night_ecbc[i] * heat_carrying_capacity_of_air / 3.1
        else:
            nightcool_energy_hourly_ecbc_value = 0
        
        nightcool_energy_hourly_ecbc.append(nightcool_energy_hourly_ecbc_value)
        nightcool_energy_ecbc = nightcool_energy_ecbc + nightcool_energy_hourly_ecbc_value
        i=i+1
    nightcool_energy_ecbc = nightcool_energy_ecbc / float(building_area)
    
    return RCDD_night_bau, RCDD_night_ecbc, nightcool_energy_bau, nightcool_energy_ecbc

  
def main():
    
    building_area=input("Enter Building area (in m2): ")
    orientation=input("Enter building orientation (E-W or N-S): ")
    building_height=input("Enter building height (in m): ")
    num_floor=input("Enter number of floors: ")
    asp_ratio=input("Enter building aspect ratio 1: ")
    WWR=input("Enter Window to Wall area Ratio (0 to 1): ")
    LPD=input("Enter Lighting Power Density (W/m2): ")
    EPD=input("Enter Equipment Power Density (W/m2): ")
    occupancy=input("Enter Building Occupancy (no. of people): ")
    evap_power = input("Enter evaporative cooler rated power: ")
    comfort_temp=input("Entre the comfort model fom the following for analysis\nASHRAE_std_55_TCM\nASHRAE_std_55_Adaptive_TCM\nIMAC_MM_NT\nIMAC_MM_UL\nIMAC_MM_LL\n: ")

    os.remove("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\results\\"+comfort_temp+".csv")
    comfort=open_file_to_write("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\results\\"+comfort_temp+".csv")
    comfort.write("city,Comfort_model,CDD as per 18.3 °C base temperature,CDD as per variable base temperature_BAU,CDD as per variable base temperature_ECBC,Residual CDD after Evaporative Pre-cooling_BAU,Residual CDD after Evaporative Pre-cooling_ECBC,Residual CDD for Evaporative Cooling_BAU,Residual CDD for Evaporative Cooling_ECBC,Residual CDD after night cooling_BAU,Residual CDD after night cooling_ECBC,Residual CDD after comfort cooling_BAU,Residual CDD after comfort cooling_ECBC,Cooling energy for fully AC_BAU,Cooling energy for fully AC_ECBC,Cooling energy after evaporative pre cooling_BAU,Cooling energy of evaporative precooler,Cooling energy after evaporative pre cooling_ECBC,Cooling energy of evaporative precooler,Cooling energy after evaporative cooling_BAU,Cooling energy of evaporative cooler_BAU,Cooling energy after evaporative cooling_ECBC,Cooling energy of evaporative cooler_ECBC,Cooling energy after night ventilation_BAU,Cooling energy after night ventilation_ECBC,Cooling energy after Comfort cooling_BAU,Cooling energy after Comfort cooling_ECBC,Lighting energy,Equipment energy\n")
    #city_name=input("ENTER CITY: ")
    for items in cities:
        
        city_name=items
        print("the city we are using is: "+str(items))
        print ("1==============================================================================================================================================================================")

        file_name=city_name+".csv"
        #print (file_name)
        f1=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\TMY Data\\ISHRAE\\CSV files_python\\"+file_name)

        DBT, dew_point, RH,atm_prs,global_horizontal_radiation,direct_normal, diffuse_horizontal, wind_speed= read_climate_variables(f1)    
        
        print ("2==============================================================================================================================================================================") 
        file_lat_log=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\TMY Data\\ISHRAE\\CSV files_python\\latitude_longitude.csv")
        latitude, longitude, roof_U_value, wall_U_value, window_U_value,window_SHGC_nonnorth, window_SHGC_north=find_latitude_longitude(file_lat_log,city_name)

        print("3===============================================================================================================================================================================")

        Twan, Twae, Twaw, Twas, Win_area_north, Win_area_east, Win_area_west, Win_area_south, Wall_area_north, Wall_area_east, Wall_area_west, Wall_area_south, Roof_area, VoP_north, VoP_east, VoP_west, VoP_south, heat_carrying_capacity_of_air, specific_heat_wall, specific_heat_roof, roof_surf_refl, wall_surf_refl, time_constant_BAU, time_constant_ECBC, thermal_capacitance = building_geometry(roof_U_value, wall_U_value, building_area, orientation, building_height, num_floor, asp_ratio, WWR)
        
        print("4===============================================================================================================================================================================")

        file_comfort=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\comfort_model.csv")
        avg_dbt=daily_average(file_comfort, DBT)

        print ("5==============================================================================================================================================================================")

        dbt_seven_day_mean=seven_day_mean(avg_dbt)
        
        print ("6==============================================================================================================================================================================")
        
        dbt_thirty_day_mean=thirty_day_mean(avg_dbt)
        
        print ("7==============================================================================================================================================================================")
        daily_mean, daily_max, daily_min = daily_mean_IMAC(DBT)
        print ("7.1==============================================================================================================================================================================")
        ASHRAE_std_55, ASHRAE_Adaptive, IMAC_MM_NT, IMAC_MM_UL, IMAC_MM_LL=comfort_model(avg_dbt, dbt_seven_day_mean, dbt_thirty_day_mean, daily_mean)
        print ("7.2==============================================================================================================================================================================")
        if comfort_temp=="ASHRAE_std_55_TCM":
            set_point = set_point_temp(ASHRAE_std_55, DBT)
        if comfort_temp=="ASHRAE_std_55_Adaptive_TCM":
            set_point = set_point_temp(ASHRAE_Adaptive, DBT)
        if comfort_temp=="IMAC_MM_NT":
            set_point = set_point_temp(IMAC_MM_NT, DBT)
        if comfort_temp=="IMAC_MM_UL":
            set_point = set_point_temp(IMAC_MM_UL, DBT)
        if comfort_temp=="IMAC_MM_LL":
            set_point = set_point_temp(IMAC_MM_LL, DBT)
       
        #set_point_temp(ASHRAE_Std_55_Adaptive_Model, IMAC_MM_Neutral_temp, IMAC_MM_Upper_limit, IMAC_MM_Lower_limit, comfort_temp)
        print ("8==============================================================================================================================================================================")
        
        file_Fabric=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\Fabric_gain.csv")	
        file_Solar=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\Solar_radiation.csv")
        
        print ("9==============================================================================================================================================================================")
        
        convert_local_Solar(longitude, file_Fabric)

        print ("10==============================================================================================================================================================================")
        file_solar_time_in_minute=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\solar_time_in_minute.txt")
        hour_angle_calculation(file_solar_time_in_minute)
        file_solar_time_in_minute.close() 
        
        print ("11==============================================================================================================================================================================")

        hour_angle_file=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\hour_angle.txt")
        do_calculations_for_solar_radiation(solar_time_in_minute, file_Solar, hour_angle_file, latitude)
        hour_angle_file.close()
        
        print ("12==============================================================================================================================================================================")

        solar_angles_file=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\solar_angles.txt")

        solar_radiation(solar_angles_file, global_horizontal_radiation, direct_normal, diffuse_horizontal, DBT)
        solar_angles_file.close()

        print ("13==============================================================================================================================================================================")
            
        #comfort_model_file=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\comfort_model.txt")
        print ("14==============================================================================================================================================================================")

        file_internal=open_file_to_read("I:\\study\\CEPT\\4th sem\\Thesis\\Arjun thesis\\Tool development\\Excel\\CSVs\\Internal_load.csv")
        occu_heat_gain, equip_heat_gain, light_heat_gain, total_internal_gain_bau, total_internal_gain_ecbc, internal_temp_rise_bau_value, internal_temp_rise_ecbc, equipment_energy, lighting_energy=internal_load(heat_carrying_capacity_of_air, file_internal, occupancy, building_area, EPD, LPD, window_SHGC_nonnorth, window_SHGC_north, Win_area_south, Win_area_west, Win_area_east, Win_area_north, solar_radiation_south, solar_radiation_east, solar_radiation_west, solar_radiation_north, solar_radiation_horizontal)
        file_internal.close()

        print ("15==============================================================================================================================================================================")

        envelope_gain(heat_carrying_capacity_of_air, building_area, building_height, set_point, Wall_area_south, Wall_area_west, Wall_area_east, Wall_area_north, Roof_area, wall_surf_refl, roof_surf_refl, wall_U_value, roof_U_value, DBT, solar_radiation_south, solar_radiation_east, solar_radiation_west, solar_radiation_north, solar_radiation_horizontal)

        print ("15.1==================================================================================================================")
        Out_humdrat, In_humdrat = humdity_ratio(DBT, RH, atm_prs, set_point)
        #print("out_humidrat"+str(out_humdrat)+"\n")

        print ("16==============================================================================================================================================================================")

        avg_night(DBT)
        
        print ("17==============================================================================================================================================================================")

        daily_avg_night(avg_night)
        
        print("18======================================================================================================================")
        
        mitigation_at_night(file_internal, time_constant_BAU, time_constant_ECBC, thermal_capacitance, heat_carrying_capacity_of_air, occupancy_schedule, set_point)

        print("19======================================================================================================================")

        fan_temp_rise=fan_gain(heat_carrying_capacity_of_air)

        print("20======================================================================================================================")

        latent_temp_rise_bau, latent_temp_rise_ecbc = latent_load(occupancy, heat_carrying_capacity_of_air, Out_humdrat, In_humdrat, fan_temp_rise, atm_prs, set_point, occupancy_schedule, total_internal_gain_bau, total_internal_gain_ecbc, internal_temp_rise_bau_value, internal_temp_rise_ecbc, mitigation_at_night_BAU, mitigation_at_night_ECBC)

        print("21======================================================================================================================")

        CDD_18_3=CDD_std(DBT)
    
        print("22======================================================================================================================")

        bpt_bau, bpt_ecbc, CDD_BAU_energy, CDD_ECBC_energy, cooling_energy_AC_bau, cooling_energy_AC_ecbc, CDD_AC_BAU, CDD_AC_ECBC= base_temperature(heat_carrying_capacity_of_air, occupancy_schedule, DBT, set_point, envelope_temperature_rise_BAU, envelope_temperature_rise_ECBC, latent_temp_rise_bau, latent_temp_rise_ecbc, mitigation_at_night_BAU, mitigation_at_night_ECBC, internal_temp_rise_bau, internal_temp_rise_ecbc, fan_temp_rise, building_area)

        print("23======================================================================================================================")

        RCDH_comfcool_bau, RCDH_comfcool_ecbc,cooling_energy_comfcool_bau, cooling_energy_comfcool_ecbc, RCDD_comfcool_bau, RCDD_comfcool_ecbc = comfort_cooling(occupancy_schedule, DBT, wind_speed, bpt_bau, bpt_ecbc, heat_carrying_capacity_of_air, building_area)

        print("24======================================================================================================================")
 
        Tw=wet_bulb(DBT,RH)

        print("24======================================================================================================================")

        pre_cool_RCDD_bau, pre_cool_RCDD_ecbc, evap_cool_RCDD_ecbc, evap_cool_RCDD_bau, evap_precool_energy_bau, evap_precool_energy_ecbc, evap_cool_energy_bau, evap_cool_energy_ecbc, evap_precooling_energy, evap_cooling_energy_bau, evap_cooling_energy_ecbc = evaporative_cooling(DBT, set_point, bpt_bau, bpt_ecbc, Tw, occupancy_schedule, heat_carrying_capacity_of_air, building_area, evap_power)

        avg_month_single, avg_jan, avg_feb, avg_mar, avg_apr, avg_may, avg_jun, avg_jul, avg_aug, avg_sep, avg_oct, avg_nov, avg_dec = monthly_average(DBT)

        avg_monthly_value = avg_monthly (avg_month_single)

        abc = monthly_average_temparature(avg_monthly_value)

        seas_avg = seasonal_average(avg_jan, avg_feb, avg_mar, avg_apr, avg_may, avg_jun, avg_jul, avg_aug, avg_sep, avg_oct, avg_nov, avg_dec, DBT)

        Min_jan, Min_feb, Min_mar, Min_apr, Min_may, Min_jun, Min_jul, Min_aug, Min_sep, Min_oct, Min_nov, Min_dec, monthly_min = monthly_minimum(DBT)

        minimum_seasonal = seasonal_minimum(Min_jan, Min_feb, Min_mar, Min_apr, Min_may, Min_jun, Min_jul, Min_aug, Min_sep, Min_oct, Min_nov, Min_dec)

        RCDD_night_bau, RCDD_night_ecbc, nightcool_energy_bau, nightcool_energy_ecbc = night_ventilation(DBT, abc, seas_avg, monthly_min, minimum_seasonal, daily_max, daily_min, avg_dbt, bpt_bau, bpt_ecbc, heat_carrying_capacity_of_air, building_area)

        
        comfort.write(str(items)+","+str(comfort_temp)+","+str(CDD_18_3)+","+str(CDD_AC_BAU)+","+str(CDD_AC_ECBC)+","+str(pre_cool_RCDD_bau)+","+str(pre_cool_RCDD_ecbc)+","+str(evap_cool_RCDD_bau)+","+str(evap_cool_RCDD_ecbc)+","+str(RCDD_night_bau)+","+str(RCDD_night_ecbc)+","+str(RCDD_comfcool_bau)+","+str(RCDD_comfcool_ecbc)+","+str(cooling_energy_AC_bau)+","+str(cooling_energy_AC_ecbc)+","+str(evap_precool_energy_bau)+","+str(evap_precooling_energy)+","+str(evap_precool_energy_ecbc)+","+str(evap_precooling_energy)+","+str(evap_cool_energy_bau)+","+str(evap_cooling_energy_bau)+","+str(evap_cool_energy_ecbc)+","+str(evap_cooling_energy_bau)+","+str(nightcool_energy_bau)+","+str(nightcool_energy_ecbc)+","+str(cooling_energy_comfcool_bau)+","+str(cooling_energy_comfcool_ecbc)+","+str(lighting_energy)+","+str(equipment_energy)+"\n")
        
        del DBT[:]
        del dew_point[:]
        del RH[:]
        del atm_prs[:]
        del global_horizontal_radiation[:]
        del direct_normal[:]
        del diffuse_horizontal[:]
        del wind_speed[:]
        del day_number[:]
        del B[:]
        del declination[:]
        del hour_angle[:]
        del date_lst[:]
        del time_lst[:]
        del avg_dbt[:]
        del occu_heat_gain[:]
        del equip_heat_gain[:]
        del light_heat_gain[:]
        del solar_time_minute[:]
        del solar_time_in_minute[:]
        del dbt_seven_day_mean[:]
        del dbt_thirty_day_mean[:]
        del daily_mean[:]
        del zenith_angle[:]
        del incident_angle_south[:]
        del incident_angle_west[:]
        del incident_angle_east[:]
        del incident_angle_north[:]
        del incident_angle_horizontal[:]
        del solar_radiation_south[:]
        del solar_radiation_west[:]
        del solar_radiation_east[:]
        del solar_radiation_north[:]
        del solar_radiation_horizontal[:]
        del dbt_temp[:]
        del ASHRAE_std_55[:]
        del ASHRAE_Adaptive[:]
        del IMAC_MM_NT[:]
        del IMAC_MM_UL[:]
        del IMAC_MM_LL[:]
        del set_point[:]
        del set_point_temperature[:]
        del radiant_gain_south_BAU[:]
        del radiant_gain_south_ECBC[:]
        del radiant_gain_west_BAU[:]
        del radiant_gain_west_ECBC[:]
        del radiant_gain_east_BAU[:]
        del radiant_gain_east_ECBC[:]
        del radiant_gain_north_BAU[:]
        del radiant_gain_north_ECBC[:]
        del radiant_gain_horizontal_BAU[:]
        del radiant_gain_horizontal_ECBC[:]
        del total_internal_gain_bau[:]
        del total_internal_gain_ecbc[:]
        del sol_air_temp_south[:]
        del sol_air_temp_west[:]
        del sol_air_temp_east[:]
        del sol_air_temp_north[:]
        del sol_air_temp_roof[:]
        del heat_gain_south_bau[:]
        del heat_gain_south_ecbc[:]
        del heat_gain_west_bau[:]
        del heat_gain_west_ecbc[:]
        del heat_gain_east_bau[:]
        del heat_gain_east_ecbc[:]
        del heat_gain_north_bau[:]
        del heat_gain_north_ecbc[:]
        del heat_gain_roof_bau[:]
        del heat_gain_roof_ecbc[:]
        del envelope_gain_bau[:]
        del envelope_gain_ecbc[:]
        del Out_humdrat[:]
        del In_humdrat[:]
        del time[:]
        del avg_night_time[:]
        del avg_night_daily[:]
        del occupancy_schedule[:]
        del lighting_schedule[:]
        del equipment_schedule[:]
        del mitigation_ECBC[:]
        del occupant_moisture_gain[:]
        del supply_air_temp_BAU[:]
        del supply_air_temp_ECBC[:]
        del supply_humdrat_bau[:]
        del supply_humdrat_ecbc[:]
        del on_coil_moisture_bau[:]
        del on_coil_moisture_ecbc[:]
        del Room_air_dew_point_bau[:]
        del Room_air_dew_point_ecbc[:]
        del infiltration[:]
        del envelope_temperature_rise_BAU[:]
        del envelope_temperature_rise_ECBC[:]
        del mitigation_at_night_BAU[:]
        del mitigation_at_night_ECBC[:]
        del latent_temp_rise_bau[:]
        del latent_temp_rise_ecbc[:]
        del internal_temp_rise_bau[:]
        del internal_temp_rise_ecbc[:]
        del bpt_bau[:]
        del bpt_ecbc[:]
        del CDH_bau[:]
        del CDH_ECBC[:]
        del RCDH_comfcool_bau[:]
        del RCDH_comfcool_ecbc[:]
        del WBT[:]
        del Tw[:]
        del RCDH_precool_bau[:]
        del RCDH_precool_ecbc[:]
        del RCDH_switch_bau[:]
        del RCDH_switch_ecbc[:]
        del avg_month_single[:]
        del avg_monthly_value[:]
        del abc[:]
        del avg_seasonal_value[:]
        del seas_avg[:]
        del min_jan[:]
        del min_feb[:]
        del min_mar[:]
        del min_apr[:]
        del min_may[:]
        del min_jun[:]
        del min_jul[:]
        del min_aug[:]
        del min_sep[:]
        del min_oct[:]
        del min_nov[:]
        del min_dec[:]
        del month_min[:]
        del monthly_min[:]
        del min_seasonal[:]
        del minimum_seasonal[:]
        del daily_max[:]
        del daily_min[:]
        del TinmaxECBC[:]
        del TinminECBC[:]
        del TinmaxBAU[:]
        del TinminBAU[:]
        del average_DBT[:]
        del daily_maximum[:]
        del daily_minimum[:]
        del amplitude_bau[:]
        del amplitude_ecbc[:]
        del sine_theeta[:]
        del angle_radian[:]
        del Tin_hourly_bau[:]
        del Tin_hourly_ecbc[:]
        del RCDH_night_bau[:]
        del RCDH_night_ecbc[:]
        del evap_precool_energy_hourly_bau[:]
        del evap_precool_energy_hourly_ecbc[:]
        del evap_cool_energy_hourly_bau[:]
        del evap_cool_energy_hourly_ecbc[:]
        del nightcool_energy_hourly_bau[:]
        del nightcool_energy_hourly_ecbc[:]
        del evap_temp[:]

if __name__ == "__main__":
    main()	
