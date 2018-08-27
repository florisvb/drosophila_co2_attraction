import analyze_shake_fly_mosquito
import figurefirst

if __name__ == '__main__':
    
    figurefirst.regenerate.clear_fifidata('supp_fig_10_shakefly_data.dillpickle', 'all')

    analyze_shake_fly_mosquito.plot_all_mosquitoes()