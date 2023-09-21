OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7282038) q[0];
sx q[0];
rz(-2.0079186) q[0];
sx q[0];
rz(1.5490305) q[0];
rz(-1.449861) q[1];
sx q[1];
rz(-2.4843042) q[1];
sx q[1];
rz(-0.54418286) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3639987) q[0];
sx q[0];
rz(-1.6382123) q[0];
sx q[0];
rz(-1.534034) q[0];
rz(-1.6149701) q[2];
sx q[2];
rz(-1.4695115) q[2];
sx q[2];
rz(-0.12461187) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6737353) q[1];
sx q[1];
rz(-2.2506672) q[1];
sx q[1];
rz(0.52635898) q[1];
rz(-pi) q[2];
rz(3.1129684) q[3];
sx q[3];
rz(-1.588436) q[3];
sx q[3];
rz(0.48354917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.071775285) q[2];
sx q[2];
rz(-1.2640307) q[2];
sx q[2];
rz(-1.7791746) q[2];
rz(-0.028256265) q[3];
sx q[3];
rz(-1.7621721) q[3];
sx q[3];
rz(2.3513667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4966999) q[0];
sx q[0];
rz(-2.4364478) q[0];
sx q[0];
rz(-1.9702966) q[0];
rz(-0.21121875) q[1];
sx q[1];
rz(-0.44208458) q[1];
sx q[1];
rz(1.7374977) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8453168) q[0];
sx q[0];
rz(-0.8527841) q[0];
sx q[0];
rz(-2.5677471) q[0];
rz(-1.4046461) q[2];
sx q[2];
rz(-1.0937905) q[2];
sx q[2];
rz(1.1039066) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6200871) q[1];
sx q[1];
rz(-1.5937935) q[1];
sx q[1];
rz(0.28316811) q[1];
rz(2.8849765) q[3];
sx q[3];
rz(-0.49391541) q[3];
sx q[3];
rz(1.4790505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.30963787) q[2];
sx q[2];
rz(-1.4761304) q[2];
sx q[2];
rz(-2.1976166) q[2];
rz(-2.5850463) q[3];
sx q[3];
rz(-0.49749938) q[3];
sx q[3];
rz(1.8054307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8495162) q[0];
sx q[0];
rz(-2.3777666) q[0];
sx q[0];
rz(-1.4235494) q[0];
rz(-2.3220093) q[1];
sx q[1];
rz(-2.3312566) q[1];
sx q[1];
rz(2.5779285) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99479988) q[0];
sx q[0];
rz(-0.99780647) q[0];
sx q[0];
rz(-1.6550199) q[0];
rz(1.7965505) q[2];
sx q[2];
rz(-0.98276897) q[2];
sx q[2];
rz(3.0622481) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3445661) q[1];
sx q[1];
rz(-2.1246506) q[1];
sx q[1];
rz(-2.7558541) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5205069) q[3];
sx q[3];
rz(-1.9120145) q[3];
sx q[3];
rz(-3.0760471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.50239572) q[2];
sx q[2];
rz(-1.1221308) q[2];
sx q[2];
rz(2.0813265) q[2];
rz(-0.075573102) q[3];
sx q[3];
rz(-2.2836756) q[3];
sx q[3];
rz(0.11463541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8183427) q[0];
sx q[0];
rz(-0.30775726) q[0];
sx q[0];
rz(2.9484205) q[0];
rz(1.5441783) q[1];
sx q[1];
rz(-1.9924106) q[1];
sx q[1];
rz(-0.65778041) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5415216) q[0];
sx q[0];
rz(-2.9642448) q[0];
sx q[0];
rz(2.970201) q[0];
x q[1];
rz(2.8305956) q[2];
sx q[2];
rz(-1.9823091) q[2];
sx q[2];
rz(-0.81931537) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.12280497) q[1];
sx q[1];
rz(-2.0969166) q[1];
sx q[1];
rz(-2.8469574) q[1];
rz(-0.25810453) q[3];
sx q[3];
rz(-1.2873642) q[3];
sx q[3];
rz(-2.5131445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0458935) q[2];
sx q[2];
rz(-0.4102439) q[2];
sx q[2];
rz(2.0945385) q[2];
rz(2.0783157) q[3];
sx q[3];
rz(-1.1626817) q[3];
sx q[3];
rz(-1.261238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3445774) q[0];
sx q[0];
rz(-2.8224967) q[0];
sx q[0];
rz(1.0239333) q[0];
rz(-2.3539885) q[1];
sx q[1];
rz(-1.5658295) q[1];
sx q[1];
rz(0.62686282) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89676566) q[0];
sx q[0];
rz(-1.5875495) q[0];
sx q[0];
rz(1.5411975) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2704029) q[2];
sx q[2];
rz(-1.8744178) q[2];
sx q[2];
rz(2.6189569) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.53351952) q[1];
sx q[1];
rz(-0.43577172) q[1];
sx q[1];
rz(2.3418952) q[1];
x q[2];
rz(0.12368006) q[3];
sx q[3];
rz(-1.9607753) q[3];
sx q[3];
rz(-2.9434413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.91288599) q[2];
sx q[2];
rz(-1.9084946) q[2];
sx q[2];
rz(-1.0181001) q[2];
rz(2.5300238) q[3];
sx q[3];
rz(-1.8194149) q[3];
sx q[3];
rz(-2.1900246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4488895) q[0];
sx q[0];
rz(-1.792181) q[0];
sx q[0];
rz(1.942379) q[0];
rz(1.9723643) q[1];
sx q[1];
rz(-1.0511845) q[1];
sx q[1];
rz(0.32454023) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51844653) q[0];
sx q[0];
rz(-1.070797) q[0];
sx q[0];
rz(2.6536233) q[0];
rz(-0.89914497) q[2];
sx q[2];
rz(-0.88485826) q[2];
sx q[2];
rz(1.9859973) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3970916) q[1];
sx q[1];
rz(-0.53452864) q[1];
sx q[1];
rz(-1.5221273) q[1];
rz(-pi) q[2];
rz(0.59431521) q[3];
sx q[3];
rz(-0.50110498) q[3];
sx q[3];
rz(-1.8584115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.67363182) q[2];
sx q[2];
rz(-1.8362074) q[2];
sx q[2];
rz(1.8661873) q[2];
rz(-1.0547137) q[3];
sx q[3];
rz(-2.349699) q[3];
sx q[3];
rz(1.5462497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4642898) q[0];
sx q[0];
rz(-1.807656) q[0];
sx q[0];
rz(1.7328847) q[0];
rz(2.6858221) q[1];
sx q[1];
rz(-2.9401638) q[1];
sx q[1];
rz(-1.2021525) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8257608) q[0];
sx q[0];
rz(-1.9681265) q[0];
sx q[0];
rz(-2.2560918) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7289594) q[2];
sx q[2];
rz(-0.79421439) q[2];
sx q[2];
rz(-2.7318294) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4533838) q[1];
sx q[1];
rz(-1.6628478) q[1];
sx q[1];
rz(-1.4870912) q[1];
x q[2];
rz(-2.8771411) q[3];
sx q[3];
rz(-1.7489986) q[3];
sx q[3];
rz(-1.4101392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.000164) q[2];
sx q[2];
rz(-1.2951853) q[2];
sx q[2];
rz(-0.84623519) q[2];
rz(2.6464461) q[3];
sx q[3];
rz(-0.64619243) q[3];
sx q[3];
rz(1.600986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071035944) q[0];
sx q[0];
rz(-1.8510011) q[0];
sx q[0];
rz(-2.0284247) q[0];
rz(-0.69560266) q[1];
sx q[1];
rz(-0.38989392) q[1];
sx q[1];
rz(1.6092469) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0198062) q[0];
sx q[0];
rz(-2.0397423) q[0];
sx q[0];
rz(-1.8437587) q[0];
x q[1];
rz(-0.74560994) q[2];
sx q[2];
rz(-0.62586212) q[2];
sx q[2];
rz(-1.9288043) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.487405) q[1];
sx q[1];
rz(-2.020917) q[1];
sx q[1];
rz(-0.022189157) q[1];
rz(-pi) q[2];
rz(1.4523776) q[3];
sx q[3];
rz(-0.97202557) q[3];
sx q[3];
rz(2.1042202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9939076) q[2];
sx q[2];
rz(-2.9116178) q[2];
sx q[2];
rz(-2.2873986) q[2];
rz(1.2285852) q[3];
sx q[3];
rz(-1.5649786) q[3];
sx q[3];
rz(1.6555697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5478741) q[0];
sx q[0];
rz(-0.49867189) q[0];
sx q[0];
rz(2.4771931) q[0];
rz(1.1876855) q[1];
sx q[1];
rz(-2.4408051) q[1];
sx q[1];
rz(-1.2845576) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5212095) q[0];
sx q[0];
rz(-1.8561583) q[0];
sx q[0];
rz(0.83831212) q[0];
rz(0.2662439) q[2];
sx q[2];
rz(-1.2028482) q[2];
sx q[2];
rz(-3.0014696) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6397275) q[1];
sx q[1];
rz(-1.1266536) q[1];
sx q[1];
rz(2.417454) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5434389) q[3];
sx q[3];
rz(-1.1809071) q[3];
sx q[3];
rz(2.011812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.61218843) q[2];
sx q[2];
rz(-0.89451423) q[2];
sx q[2];
rz(-0.55076304) q[2];
rz(2.4380056) q[3];
sx q[3];
rz(-1.0102605) q[3];
sx q[3];
rz(1.5065058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4187014) q[0];
sx q[0];
rz(-2.7846865) q[0];
sx q[0];
rz(-0.044145949) q[0];
rz(1.528953) q[1];
sx q[1];
rz(-1.9020558) q[1];
sx q[1];
rz(2.3619161) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75694376) q[0];
sx q[0];
rz(-1.0352408) q[0];
sx q[0];
rz(1.4063565) q[0];
x q[1];
rz(-0.463562) q[2];
sx q[2];
rz(-0.87052411) q[2];
sx q[2];
rz(-2.314032) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.838678) q[1];
sx q[1];
rz(-2.089114) q[1];
sx q[1];
rz(1.360449) q[1];
rz(-pi) q[2];
x q[2];
rz(2.036282) q[3];
sx q[3];
rz(-2.2959384) q[3];
sx q[3];
rz(-2.5940965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.49446517) q[2];
sx q[2];
rz(-0.72138849) q[2];
sx q[2];
rz(-1.54281) q[2];
rz(0.70458448) q[3];
sx q[3];
rz(-1.6274118) q[3];
sx q[3];
rz(3.1295479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4298532) q[0];
sx q[0];
rz(-1.5562417) q[0];
sx q[0];
rz(-0.65162311) q[0];
rz(-1.7977057) q[1];
sx q[1];
rz(-1.6603036) q[1];
sx q[1];
rz(2.4659326) q[1];
rz(-0.96991878) q[2];
sx q[2];
rz(-2.6261332) q[2];
sx q[2];
rz(0.020608227) q[2];
rz(-2.5034954) q[3];
sx q[3];
rz(-1.4641855) q[3];
sx q[3];
rz(-2.2104213) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];