import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

class Platte:
    def __init__(self, filename_innerforce, filename_raction, d_mm, c_nom_mm, d4_mm, a4_mm, d3_mm, a3_mm, s_mm, f_cd_Nmm2, tau_cd_Nmm2, Dmax_mm, f_sd_Nmm2, E_s_Nmm2, Richtung):
        self.Plattendicke_d_mm = d_mm
        self.Ueberdeckung_cnom_mm = c_nom_mm
        self.Vierte_Lage_d_mm, self.Teilung_Vierte_Lage_a_mm = d4_mm, a4_mm
        self.Dritte_Lage_d_mm, self.Teilung_Dritte_Lage_a_mm = d3_mm, a3_mm
        self.Einsenkung_s_mm = s_mm
        self.Druckfestigkeit_f_cd_Nmm2 = f_cd_Nmm2
        self.Schubfestigkeit_tau_cd_Nmm2 = tau_cd_Nmm2
        self.Maximalkorn_dmax_mm = Dmax_mm
        self.Stahlfliessgrenze_fsd_Nmm2 = f_sd_Nmm2
        self.Stahlsteifigkeit_E_s_Nmm2 = E_s_Nmm2
        self.Ausrichtung_4_Lage = Richtung

        self.inner_force = pd.read_csv(filename_innerforce, names=[
            'Richt. x', 'mx[kN]', 'my[kN]', 'mxy[kN]', 'vx[kN/m]', 'vy[kN/m]', 'm1[kN]', 'm2[kN]', 'Winkel[°]'
            ], sep=';').dropna(axis=0).reset_index()
        
        self.reaction = pd.read_csv(filename_raction, names=[
            'x[m]', 'y[m]', 'rZmin[kN/m2]', 'rZmax[kN/m2]'
            ], sep=';')
        
        self.connected_data = pd.concat([self.reaction[['x[m]', 'y[m]']], self.inner_force[['vx[kN/m]', 'vy[kN/m]', 'mx[kN]', 'my[kN]', 'mxy[kN]']]], axis=1)
        
        self.statische_Hoehe_d3_mm = self.Plattendicke_d_mm - self.Ueberdeckung_cnom_mm - self.Vierte_Lage_d_mm - self.Dritte_Lage_d_mm/2
        self.statische_Hoehe_d4_mm = self.Plattendicke_d_mm - self.Ueberdeckung_cnom_mm - self.Vierte_Lage_d_mm/2
        self.statische_Hoehe_dv_mm = (self.statische_Hoehe_d4_mm + self.statische_Hoehe_d3_mm) / 2 - self.Einsenkung_s_mm


    def get_filename(self, Typ):
        filename = Typ + '_' + 'Platte_d_' + str(self.Plattendicke_d_mm) + '_c_nom_' + str(
        self.Ueberdeckung_cnom_mm) + '_1_Lage_' + self.Ausrichtung_4_Lage + '_d4_' + str(
        self.Vierte_Lage_d_mm) + '_' + str(
        self.Teilung_Vierte_Lage_a_mm) + '_d3_' + str(self.Dritte_Lage_d_mm) + '_' + str(
        self.Teilung_Dritte_Lage_a_mm) + '_' + self.str_date()
        return filename
    
    def save_csv(self, dataframe, Typ):
        filename = self.get_filename(Typ=Typ)
        dataframe.to_csv(filename + '.csv')

    def str_date(self):
        now = datetime.now() # current date and time
        year = now.strftime("%Y")
        month = now.strftime("%m")
        day = now.strftime("%d")
        hour = now.strftime("%H")
        minute = now.strftime("%M")
        second = now.strftime("%S")
        date = year + '_' + month + '_' + day + '__' + hour + '_' + minute + '_' + second
        return date

    def m_4_Rd(self):
        m_vier_Rd = -(self.Vierte_Lage_d_mm**2*np.pi/4*1000/self.Teilung_Vierte_Lage_a_mm)*self.Stahlfliessgrenze_fsd_Nmm2*(
            (self.Plattendicke_d_mm-self.Ueberdeckung_cnom_mm-self.Vierte_Lage_d_mm/2-self.Einsenkung_s_mm)-(
            self.Vierte_Lage_d_mm**2*np.pi/4*1000/self.Teilung_Vierte_Lage_a_mm)*self.Stahlfliessgrenze_fsd_Nmm2/(
            2*1000*self.Druckfestigkeit_f_cd_Nmm2))/10**6
        return m_vier_Rd

    def m_3_Rd(self):
        m_drei_Rd = -(self.Dritte_Lage_d_mm**2*np.pi/4*1000/self.Teilung_Dritte_Lage_a_mm)*self.Stahlfliessgrenze_fsd_Nmm2*(
            (self.Plattendicke_d_mm-self.Ueberdeckung_cnom_mm-self.Vierte_Lage_d_mm-self.Dritte_Lage_d_mm/2-self.Einsenkung_s_mm)-(
            self.Dritte_Lage_d_mm**2*np.pi/4*1000/self.Teilung_Dritte_Lage_a_mm)*self.Stahlfliessgrenze_fsd_Nmm2/(
            2*1000*self.Druckfestigkeit_f_cd_Nmm2))/10**6
        return m_drei_Rd
    
    def m_x_Rd(self):
        if self.Ausrichtung_4_Lage == 'x':
            return self.m_4_Rd()
        else:
            return self.m_3_Rd()
        
    def m_y_Rd(self):
        if self.Ausrichtung_4_Lage == 'x':
            return self.m_3_Rd()
        else:
            return self.m_4_Rd()
                   
    def Biegenachweise(self):
        result = self.connected_data.loc[:,['x[m]', 'y[m]', 'vx[kN/m]', 'vy[kN/m]', 'mx[kN]', 'my[kN]', 'mxy[kN]']]
        result['phi_0[rad]'] = np.arctan(self.inner_force['vy[kN/m]']/self.inner_force['vx[kN/m]'])
        result['m_n_d[kN]'] = self.inner_force['mx[kN]']*np.cos(result['phi_0[rad]'])**2 + self.inner_force['my[kN]']*np.sin(
            result['phi_0[rad]'])**2 + self.inner_force['mxy[kN]']*np.sin(2*result['phi_0[rad]'])
        result['m_n_Rd[kN]'] = self.m_x_Rd() * np.cos(result['phi_0[rad]'])**2 + self.m_y_Rd() * np.sin(result['phi_0[rad]'])**2
        result['m_n_d/m_n_Rd[-]'] = abs(result['m_n_d[kN]'] / result['m_n_Rd[kN]'])
        return result

    def Querkraftnachweise(self):
        result = self.Biegenachweise()
        result['zeta[-]'] = 1 / (np.sin(result['phi_0[rad]'])**4 + np.cos(result['phi_0[rad]'])**4)
        result['e_v[10^-3]'] = self.Stahlfliessgrenze_fsd_Nmm2 / self.Stahlsteifigkeit_E_s_Nmm2 * result['m_n_d/m_n_Rd[-]'] / 10**(-3)
        result['k_g[-]'] = 48 / (16 + self.Maximalkorn_dmax_mm)
        result['k_d[-]'] = 1 / (1 + result['e_v[10^-3]']*10**(-3)*(self.statische_Hoehe_dv_mm + self.Einsenkung_s_mm) * result['k_g[-]'] * result['zeta[-]'])
        result['v0d[kN/m]'] = (result['vx[kN/m]']**2 + result['vy[kN/m]']**2)**(1/2)
        result['v_Rd[kN/m]'] = result['k_d[-]'] * self.Schubfestigkeit_tau_cd_Nmm2 * self.statische_Hoehe_dv_mm
        result['v_d / v_Rd[-]'] = result['v0d[kN/m]'] / result['v_Rd[kN/m]']
        return result
       
    def Querkraftnachweise_nicht_erfuellt(self, sort_values=True, save_csv=False):
        result = self.Querkraftnachweise()
        V_bad = result.loc[result['v_d / v_Rd[-]'] > 1.00]
        if sort_values:
            V_bad = V_bad.sort_values(by='v_d / v_Rd[-]', ascending=False)
        if save_csv:
            self.save_csv(dataframe=V_bad, Typ='V')
        return V_bad
    
    def Biegenachweise_nicht_erfuellt(self, sort_values=True, save_csv=False):
        result = self.Biegenachweise()
        M_bad = result.loc[result['m_n_d/m_n_Rd[-]'] > 1.00]
        if sort_values:
             M_bad = M_bad.sort_values(by='m_n_d/m_n_Rd[-]', ascending=False)
        if save_csv:
            self.save_csv(dataframe=M_bad, Typ='M')
        return M_bad
    
    def Querkraftnachweise_erfuellt(self, sort_values=True, save_csv=False):
        result = self.Querkraftnachweise()
        V_good = result.loc[result['v_d / v_Rd[-]'] <= 1.00]
        if sort_values:
            V_good = V_good.sort_values(by='v_d / v_Rd[-]', ascending=False)
        if save_csv:
            self.save_csv(dataframe=V_good, Typ='V')
        return V_good
    
    def Biegenachweise_erfuellt(self, sort_values=True, save_csv=False):
        result = self.Biegenachweise()
        M_good = result.loc[result['m_n_d/m_n_Rd[-]'] <= 1.00]
        if sort_values:
             M_good = M_good.sort_values(by='m_n_d/m_n_Rd[-]', ascending=False)
        if save_csv:
            self.save_csv(dataframe=M_good, Typ='V')
        return M_good
    
    def print_info(self):
        if self.Ausrichtung_4_Lage == 'x':
            d_x, a_x = self.Vierte_Lage_d_mm, self.Teilung_Vierte_Lage_a_mm
            d_y, a_y = self.Dritte_Lage_d_mm, self.Teilung_Dritte_Lage_a_mm
            d_x_s = self.statische_Hoehe_d4_mm 
            d_y_s = self.statische_Hoehe_d3_mm 
        else:
            d_x, a_x = self.Dritte_Lage_d_mm, self.Teilung_Dritte_Lage_a_mm
            d_y, a_y = self.Vierte_Lage_d_mm, self.Teilung_Vierte_Lage_a_mm
            d_x_s = self.statische_Hoehe_d3_mm 
            d_y_s = self.statische_Hoehe_d4_mm 

        def V_Rd_min():
            zeta = 2.0
            tau_cd = self.Schubfestigkeit_tau_cd_Nmm2
            e_v = self.Stahlfliessgrenze_fsd_Nmm2 / self.Stahlsteifigkeit_E_s_Nmm2
            k_g = 48 / (16 + self.Maximalkorn_dmax_mm)
            d_v = self.statische_Hoehe_dv_mm
            d = self.statische_Hoehe_d4_mm
            k_d = 1 / (1 + e_v * d * k_g * zeta)
            V_Rd_min = k_d * tau_cd * d_v
            return round(V_Rd_min, 2)
        
        M_failed = round(len(self.Biegenachweise_nicht_erfuellt()) / len(self.Biegenachweise()) * 100, 2)
        V_failed = round(len(self.Querkraftnachweise_nicht_erfuellt()) / len(self.Querkraftnachweise()) * 100, 2)

        print('Decke d =', str(self.Plattendicke_d_mm) + ' mm, c_nom = ' + str(self.Ueberdeckung_cnom_mm) 
              + ' mm, Einsenkung s = ' + str(self.Einsenkung_s_mm) + ' mm')
        print('---------------------------------------------------------------------------------------------------')
        print('Bewehrung: x ' + str(d_x) + '|' + str(a_x) + ' mm, y ' + str(d_y) + '|' 
              + str(a_y) + ' mm' + ', f_sd = ' + str(self.Stahlfliessgrenze_fsd_Nmm2) + ' N/mm2, E_s = ' + str(self.Stahlsteifigkeit_E_s_Nmm2) + ' N/mm2')
        print('Beton: f_cd = ' + str(self.Druckfestigkeit_f_cd_Nmm2) + ' N/mm2, tau_cd = ' + str(self.Schubfestigkeit_tau_cd_Nmm2) + '0 N/mm2, D_max = ' + str(self.Maximalkorn_dmax_mm) + ' mm')
        print('---------------------------------------------------------------------------------------------------')
        print('m_x_Rd(d_x = ' + str(d_x_s) + ' mm) = ' +  str(round(self.m_x_Rd(), 2)) + ' kNm/m')
        print('m_y_Rd(d_y = ' + str(d_y_s) + ' mm) = ' +  str(round(self.m_y_Rd(), 2)) + ' kNm/m')
        print('v_Rd_min = ' + str(V_Rd_min()) + ' kN/m (für Nachweise unten irrelevant)')
        print('d_v = ' + str(self.statische_Hoehe_dv_mm) + ' mm')
        print('---------------------------------------------------------------------------------------------------')
        print('Gescheiterte Biegenachweise: ' + str(len(self.Biegenachweise_nicht_erfuellt())) + '/' + str(len(self.Biegenachweise())) + ' (' + str(M_failed) + ' %)')
        print('Gescheiterte Querkraftnachweise: ' + str(len(self.Querkraftnachweise_nicht_erfuellt())) + '/' + str(len(self.Querkraftnachweise()))+ ' (' + str(V_failed) + ' %)')


    def print_Nachweise(self, Typ='V', save_file=False, markersize=5, bad_only=False, annotation=False, annotation_size=5):

        if Typ == 'V':
            result_column = 'v_d / v_Rd[-]'
            latex_name = 'v_d / v_{Rd}'
            nicht_erfuellt = self.Querkraftnachweise_nicht_erfuellt()
            erfuellt = self.Querkraftnachweise_erfuellt()
        elif Typ == 'M':
            result_column = 'm_n_d/m_n_Rd[-]'
            latex_name = 'm_{nd} / m_{nRd}'
            nicht_erfuellt = self.Biegenachweise_nicht_erfuellt()
            erfuellt = self.Biegenachweise_erfuellt()

        fig, ax = plt.subplots()

        x_nicht_erfuellt = nicht_erfuellt['x[m]']
        y_nicht_erfuellt = nicht_erfuellt['y[m]']

        x_erfuellt = erfuellt['x[m]']
        y_erfuellt = erfuellt['y[m]']

        # Use the same colormap as the scatter points
        cmap = plt.cm.get_cmap('RdYlGn')
        
        if bad_only == False:
            ax.scatter(x=x_erfuellt, y=y_erfuellt, s=markersize, c='gainsboro', marker='s')

        # Update the annotation color based on the result column values
        annotations_color = cmap(1/nicht_erfuellt[result_column])

        if annotation == False:
            ax.scatter(x=x_nicht_erfuellt, y=y_nicht_erfuellt, s=markersize, c=1/nicht_erfuellt[result_column], cmap='RdYlGn', marker='s')

        if annotation:
            # Annotation
            for index in range(len(nicht_erfuellt)):
                text, x, y = round(nicht_erfuellt.iloc[index, -1], 2), nicht_erfuellt.iloc[index, 0], nicht_erfuellt.iloc[index, 1]
                color = annotations_color[index]  # Use the color from the colormap
                ax.annotate(text=str(text), xy=[x, y], fontsize=annotation_size, ha='center', c=color)
        
        ax.set_aspect('equal')

        if len(nicht_erfuellt) > 0:
            str_title = 'Nicht erfüllte Nachweise: ' + str(
                round(nicht_erfuellt[result_column].min(), 2)) + ' < $ ' + latex_name + '\leq $' + str(
                round(nicht_erfuellt[result_column].max(), 2))
        else:
            str_title = 'Alle Nachweise erfüllt: ' + str(
                round(erfuellt[result_column].min(), 2)) + ' < $ ' + latex_name + '\leq $' + str(
                round(erfuellt[result_column].max(), 2))

        ax.set_title(str_title)

        if save_file:
            filename = self.get_filename(Typ=Typ)
            plt.savefig(format='png', fname=filename + '.png', dpi=600)
            plt.savefig(format='pdf', fname=filename + '.pdf')

        plt.show()