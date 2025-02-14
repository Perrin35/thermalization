OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.74855411) q[0];
sx q[0];
rz(-1.3286123) q[0];
sx q[0];
rz(0.42088977) q[0];
rz(2.4731877) q[1];
sx q[1];
rz(-1.547812) q[1];
sx q[1];
rz(1.1512383) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6164857) q[0];
sx q[0];
rz(-0.49607222) q[0];
sx q[0];
rz(0.43323364) q[0];
rz(-2.6815979) q[2];
sx q[2];
rz(-1.056571) q[2];
sx q[2];
rz(-2.0005039) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9734263) q[1];
sx q[1];
rz(-0.81422537) q[1];
sx q[1];
rz(2.3553215) q[1];
rz(-pi) q[2];
rz(2.8731396) q[3];
sx q[3];
rz(-2.2055948) q[3];
sx q[3];
rz(-2.6688547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.067817299) q[2];
sx q[2];
rz(-2.1442118) q[2];
sx q[2];
rz(2.1201521) q[2];
rz(-1.8197618) q[3];
sx q[3];
rz(-0.50522155) q[3];
sx q[3];
rz(0.57636133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3385056) q[0];
sx q[0];
rz(-2.1144688) q[0];
sx q[0];
rz(-2.8402253) q[0];
rz(-1.1400247) q[1];
sx q[1];
rz(-2.3224484) q[1];
sx q[1];
rz(1.0637306) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50142215) q[0];
sx q[0];
rz(-2.8482264) q[0];
sx q[0];
rz(1.1430995) q[0];
rz(-pi) q[1];
rz(0.76164328) q[2];
sx q[2];
rz(-1.7178665) q[2];
sx q[2];
rz(2.0936793) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1769907) q[1];
sx q[1];
rz(-2.1637205) q[1];
sx q[1];
rz(-1.9511392) q[1];
x q[2];
rz(-2.1065936) q[3];
sx q[3];
rz(-1.4685653) q[3];
sx q[3];
rz(-0.4680948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2782044) q[2];
sx q[2];
rz(-0.74625838) q[2];
sx q[2];
rz(-2.9928652) q[2];
rz(-1.1778098) q[3];
sx q[3];
rz(-1.9494467) q[3];
sx q[3];
rz(1.8671794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3026368) q[0];
sx q[0];
rz(-1.1936854) q[0];
sx q[0];
rz(1.0021915) q[0];
rz(-1.8110555) q[1];
sx q[1];
rz(-1.9621153) q[1];
sx q[1];
rz(-0.60752121) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1992545) q[0];
sx q[0];
rz(-0.94616854) q[0];
sx q[0];
rz(-2.6074431) q[0];
rz(1.9776865) q[2];
sx q[2];
rz(-1.3156789) q[2];
sx q[2];
rz(1.6734378) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3866736) q[1];
sx q[1];
rz(-2.6133839) q[1];
sx q[1];
rz(-2.9186503) q[1];
rz(-2.1488299) q[3];
sx q[3];
rz(-0.93621636) q[3];
sx q[3];
rz(-1.6333333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9017631) q[2];
sx q[2];
rz(-0.87368691) q[2];
sx q[2];
rz(-2.152781) q[2];
rz(-1.4766988) q[3];
sx q[3];
rz(-1.3380545) q[3];
sx q[3];
rz(-0.57507676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7779509) q[0];
sx q[0];
rz(-2.2626484) q[0];
sx q[0];
rz(-1.2203891) q[0];
rz(-1.8266504) q[1];
sx q[1];
rz(-1.5053791) q[1];
sx q[1];
rz(3.0016518) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51516384) q[0];
sx q[0];
rz(-2.4068711) q[0];
sx q[0];
rz(0.61798851) q[0];
rz(-pi) q[1];
rz(2.1084505) q[2];
sx q[2];
rz(-2.2058479) q[2];
sx q[2];
rz(0.5039353) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8144758) q[1];
sx q[1];
rz(-0.45683858) q[1];
sx q[1];
rz(1.4903699) q[1];
rz(-pi) q[2];
x q[2];
rz(0.071752944) q[3];
sx q[3];
rz(-1.9535616) q[3];
sx q[3];
rz(1.07384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4466897) q[2];
sx q[2];
rz(-2.0835154) q[2];
sx q[2];
rz(-3.0042082) q[2];
rz(-1.8047699) q[3];
sx q[3];
rz(-1.5788014) q[3];
sx q[3];
rz(-0.00060753879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7430275) q[0];
sx q[0];
rz(-1.407857) q[0];
sx q[0];
rz(-2.9803357) q[0];
rz(-2.297961) q[1];
sx q[1];
rz(-0.7470986) q[1];
sx q[1];
rz(1.3074494) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3056934) q[0];
sx q[0];
rz(-0.81863716) q[0];
sx q[0];
rz(-2.7710415) q[0];
rz(-pi) q[1];
rz(2.0548115) q[2];
sx q[2];
rz(-0.40905158) q[2];
sx q[2];
rz(-1.7526) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2469663) q[1];
sx q[1];
rz(-1.481253) q[1];
sx q[1];
rz(0.98693165) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90353538) q[3];
sx q[3];
rz(-2.1318815) q[3];
sx q[3];
rz(2.4172731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1696986) q[2];
sx q[2];
rz(-2.7104388) q[2];
sx q[2];
rz(2.2141854) q[2];
rz(-1.5291322) q[3];
sx q[3];
rz(-1.1011139) q[3];
sx q[3];
rz(-2.7839938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9277495) q[0];
sx q[0];
rz(-2.1723211) q[0];
sx q[0];
rz(1.4056322) q[0];
rz(-1.4145781) q[1];
sx q[1];
rz(-2.3443293) q[1];
sx q[1];
rz(-2.3419211) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41922955) q[0];
sx q[0];
rz(-1.5088313) q[0];
sx q[0];
rz(3.1298766) q[0];
rz(-1.437101) q[2];
sx q[2];
rz(-1.8779199) q[2];
sx q[2];
rz(1.6973059) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1649014) q[1];
sx q[1];
rz(-1.3200687) q[1];
sx q[1];
rz(-1.9216955) q[1];
x q[2];
rz(0.42678435) q[3];
sx q[3];
rz(-2.7023661) q[3];
sx q[3];
rz(-0.78204417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.25524011) q[2];
sx q[2];
rz(-2.4031874) q[2];
sx q[2];
rz(-0.8026455) q[2];
rz(0.32043996) q[3];
sx q[3];
rz(-1.8010062) q[3];
sx q[3];
rz(1.6360487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36444148) q[0];
sx q[0];
rz(-0.027712263) q[0];
sx q[0];
rz(0.030315422) q[0];
rz(-1.3069356) q[1];
sx q[1];
rz(-1.0825284) q[1];
sx q[1];
rz(-2.511715) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4277074) q[0];
sx q[0];
rz(-2.6763335) q[0];
sx q[0];
rz(-0.42963877) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81923072) q[2];
sx q[2];
rz(-0.60556817) q[2];
sx q[2];
rz(1.8964963) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4640208) q[1];
sx q[1];
rz(-1.3673165) q[1];
sx q[1];
rz(-2.8427441) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9050352) q[3];
sx q[3];
rz(-0.45522296) q[3];
sx q[3];
rz(-0.61456028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.031875413) q[2];
sx q[2];
rz(-2.4745291) q[2];
sx q[2];
rz(1.9943705) q[2];
rz(-2.4540497) q[3];
sx q[3];
rz(-2.5496428) q[3];
sx q[3];
rz(-1.3233002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3654093) q[0];
sx q[0];
rz(-2.9054346) q[0];
sx q[0];
rz(2.6614406) q[0];
rz(2.7393553) q[1];
sx q[1];
rz(-1.6419342) q[1];
sx q[1];
rz(-0.38280907) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7744336) q[0];
sx q[0];
rz(-1.2437789) q[0];
sx q[0];
rz(2.9190382) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8620969) q[2];
sx q[2];
rz(-1.6812107) q[2];
sx q[2];
rz(-1.8556653) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4856163) q[1];
sx q[1];
rz(-1.1895863) q[1];
sx q[1];
rz(-0.31869546) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6005101) q[3];
sx q[3];
rz(-1.7669356) q[3];
sx q[3];
rz(-2.8063584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7499034) q[2];
sx q[2];
rz(-2.4112371) q[2];
sx q[2];
rz(-2.7833617) q[2];
rz(-2.7273438) q[3];
sx q[3];
rz(-1.0284938) q[3];
sx q[3];
rz(-2.7481368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2888032) q[0];
sx q[0];
rz(-1.2816592) q[0];
sx q[0];
rz(-2.7237256) q[0];
rz(1.6425543) q[1];
sx q[1];
rz(-2.7211029) q[1];
sx q[1];
rz(2.5299759) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8369524) q[0];
sx q[0];
rz(-1.3654183) q[0];
sx q[0];
rz(-2.9083328) q[0];
x q[1];
rz(0.29199227) q[2];
sx q[2];
rz(-0.21459178) q[2];
sx q[2];
rz(2.8923182) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4222718) q[1];
sx q[1];
rz(-0.55972717) q[1];
sx q[1];
rz(2.1156963) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2232333) q[3];
sx q[3];
rz(-2.7866552) q[3];
sx q[3];
rz(2.8508773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.186782) q[2];
sx q[2];
rz(-1.6839226) q[2];
sx q[2];
rz(-0.35624722) q[2];
rz(-0.48480836) q[3];
sx q[3];
rz(-1.3508513) q[3];
sx q[3];
rz(-3.1118605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32271069) q[0];
sx q[0];
rz(-1.1047381) q[0];
sx q[0];
rz(-1.6792962) q[0];
rz(2.7541584) q[1];
sx q[1];
rz(-2.2144364) q[1];
sx q[1];
rz(-2.3364054) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1772386) q[0];
sx q[0];
rz(-2.2573009) q[0];
sx q[0];
rz(-2.3330188) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3024955) q[2];
sx q[2];
rz(-0.57575127) q[2];
sx q[2];
rz(-2.8481399) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2428875) q[1];
sx q[1];
rz(-1.4243898) q[1];
sx q[1];
rz(0.28621788) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58918192) q[3];
sx q[3];
rz(-2.7134827) q[3];
sx q[3];
rz(-2.5758139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.33599535) q[2];
sx q[2];
rz(-0.62312859) q[2];
sx q[2];
rz(-1.9192637) q[2];
rz(-0.81021106) q[3];
sx q[3];
rz(-2.2170292) q[3];
sx q[3];
rz(1.5789825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41351086) q[0];
sx q[0];
rz(-1.6376729) q[0];
sx q[0];
rz(-1.6117657) q[0];
rz(1.275508) q[1];
sx q[1];
rz(-2.7559912) q[1];
sx q[1];
rz(-1.4082946) q[1];
rz(-1.1861943) q[2];
sx q[2];
rz(-1.6976962) q[2];
sx q[2];
rz(-3.1023956) q[2];
rz(1.0216711) q[3];
sx q[3];
rz(-2.3952978) q[3];
sx q[3];
rz(1.8412347) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
