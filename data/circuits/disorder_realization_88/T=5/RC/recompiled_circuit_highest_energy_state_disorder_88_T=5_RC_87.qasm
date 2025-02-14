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
rz(-2.5992844) q[0];
sx q[0];
rz(-3.0071654) q[0];
sx q[0];
rz(-2.0943213) q[0];
rz(-0.32416999) q[1];
sx q[1];
rz(2.9919762) q[1];
sx q[1];
rz(11.621578) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0271929) q[0];
sx q[0];
rz(-1.2978221) q[0];
sx q[0];
rz(1.2370212) q[0];
rz(-pi) q[1];
rz(-0.23350291) q[2];
sx q[2];
rz(-2.2036607) q[2];
sx q[2];
rz(-1.3695182) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8103247) q[1];
sx q[1];
rz(-2.4400418) q[1];
sx q[1];
rz(0.74358209) q[1];
rz(1.3421435) q[3];
sx q[3];
rz(-0.49352577) q[3];
sx q[3];
rz(2.760104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2904539) q[2];
sx q[2];
rz(-1.6893427) q[2];
sx q[2];
rz(-1.5645082) q[2];
rz(2.5751298) q[3];
sx q[3];
rz(-2.5201859) q[3];
sx q[3];
rz(1.2982781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8241149) q[0];
sx q[0];
rz(-0.7190187) q[0];
sx q[0];
rz(-2.4216006) q[0];
rz(-2.8672245) q[1];
sx q[1];
rz(-0.73148483) q[1];
sx q[1];
rz(-0.55363399) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0087411441) q[0];
sx q[0];
rz(-2.7692502) q[0];
sx q[0];
rz(-1.2471863) q[0];
x q[1];
rz(2.4963412) q[2];
sx q[2];
rz(-0.49075365) q[2];
sx q[2];
rz(3.0440999) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.95018847) q[1];
sx q[1];
rz(-2.2850415) q[1];
sx q[1];
rz(-0.49226239) q[1];
x q[2];
rz(1.0721562) q[3];
sx q[3];
rz(-2.0387531) q[3];
sx q[3];
rz(-1.4730723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9686034) q[2];
sx q[2];
rz(-1.9393238) q[2];
sx q[2];
rz(-2.7776862) q[2];
rz(3.0299752) q[3];
sx q[3];
rz(-2.2026187) q[3];
sx q[3];
rz(-2.3811471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97742057) q[0];
sx q[0];
rz(-0.82472473) q[0];
sx q[0];
rz(-0.0053996276) q[0];
rz(0.81575704) q[1];
sx q[1];
rz(-1.1382256) q[1];
sx q[1];
rz(0.38591787) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49026981) q[0];
sx q[0];
rz(-1.2322591) q[0];
sx q[0];
rz(-1.2744034) q[0];
rz(2.1109525) q[2];
sx q[2];
rz(-2.5084426) q[2];
sx q[2];
rz(-1.4667222) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2066127) q[1];
sx q[1];
rz(-2.7394208) q[1];
sx q[1];
rz(1.241317) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1564527) q[3];
sx q[3];
rz(-1.4477535) q[3];
sx q[3];
rz(0.70085722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.58891121) q[2];
sx q[2];
rz(-0.83325714) q[2];
sx q[2];
rz(-1.7595278) q[2];
rz(0.16119334) q[3];
sx q[3];
rz(-1.4313982) q[3];
sx q[3];
rz(0.62895697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2564119) q[0];
sx q[0];
rz(-1.2762524) q[0];
sx q[0];
rz(1.8026344) q[0];
rz(1.3171875) q[1];
sx q[1];
rz(-2.1664797) q[1];
sx q[1];
rz(0.29979527) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9343963) q[0];
sx q[0];
rz(-0.75524932) q[0];
sx q[0];
rz(-0.51934262) q[0];
x q[1];
rz(-2.1607397) q[2];
sx q[2];
rz(-1.9741657) q[2];
sx q[2];
rz(-1.9669176) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0865752) q[1];
sx q[1];
rz(-0.37913943) q[1];
sx q[1];
rz(-1.9490446) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13238975) q[3];
sx q[3];
rz(-2.2023938) q[3];
sx q[3];
rz(-1.6211444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.51912159) q[2];
sx q[2];
rz(-2.2032479) q[2];
sx q[2];
rz(2.8224714) q[2];
rz(-0.090911344) q[3];
sx q[3];
rz(-2.9967283) q[3];
sx q[3];
rz(-1.3566141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7262481) q[0];
sx q[0];
rz(-2.321796) q[0];
sx q[0];
rz(0.0023284624) q[0];
rz(-1.5709411) q[1];
sx q[1];
rz(-0.80051533) q[1];
sx q[1];
rz(-1.9470107) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9228421) q[0];
sx q[0];
rz(-0.19521579) q[0];
sx q[0];
rz(-2.4890635) q[0];
x q[1];
rz(1.5244687) q[2];
sx q[2];
rz(-2.4232658) q[2];
sx q[2];
rz(-2.3026932) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3227685) q[1];
sx q[1];
rz(-1.4684825) q[1];
sx q[1];
rz(-1.9167856) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8525328) q[3];
sx q[3];
rz(-2.1093371) q[3];
sx q[3];
rz(-2.0592214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2583367) q[2];
sx q[2];
rz(-2.363435) q[2];
sx q[2];
rz(-0.077433132) q[2];
rz(0.94610131) q[3];
sx q[3];
rz(-0.48282048) q[3];
sx q[3];
rz(2.6879123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2133863) q[0];
sx q[0];
rz(-0.086406924) q[0];
sx q[0];
rz(-1.4148096) q[0];
rz(2.4746223) q[1];
sx q[1];
rz(-1.357115) q[1];
sx q[1];
rz(-2.736843) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20807438) q[0];
sx q[0];
rz(-0.41941038) q[0];
sx q[0];
rz(2.389877) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26164345) q[2];
sx q[2];
rz(-0.930013) q[2];
sx q[2];
rz(-1.4705758) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7915899) q[1];
sx q[1];
rz(-1.5096997) q[1];
sx q[1];
rz(-1.1076124) q[1];
rz(-pi) q[2];
rz(0.36207288) q[3];
sx q[3];
rz(-2.0303465) q[3];
sx q[3];
rz(-1.5611695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3177967) q[2];
sx q[2];
rz(-2.3522289) q[2];
sx q[2];
rz(-2.6981567) q[2];
rz(-0.41380841) q[3];
sx q[3];
rz(-0.2897073) q[3];
sx q[3];
rz(-2.2558291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39171788) q[0];
sx q[0];
rz(-1.7789919) q[0];
sx q[0];
rz(0.66746563) q[0];
rz(-3.093847) q[1];
sx q[1];
rz(-2.0011438) q[1];
sx q[1];
rz(-2.023229) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55505468) q[0];
sx q[0];
rz(-1.0411052) q[0];
sx q[0];
rz(2.0425969) q[0];
x q[1];
rz(1.3862382) q[2];
sx q[2];
rz(-0.58635752) q[2];
sx q[2];
rz(1.9826012) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.075526) q[1];
sx q[1];
rz(-2.3956163) q[1];
sx q[1];
rz(-0.64167185) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0736865) q[3];
sx q[3];
rz(-2.0693681) q[3];
sx q[3];
rz(-0.25085051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.095801) q[2];
sx q[2];
rz(-1.2597151) q[2];
sx q[2];
rz(-0.0087139159) q[2];
rz(-1.2039315) q[3];
sx q[3];
rz(-0.34398505) q[3];
sx q[3];
rz(3.1194527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0225723) q[0];
sx q[0];
rz(-2.75596) q[0];
sx q[0];
rz(2.2517396) q[0];
rz(-1.285137) q[1];
sx q[1];
rz(-1.0533918) q[1];
sx q[1];
rz(3.0056675) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0014627731) q[0];
sx q[0];
rz(-0.23389947) q[0];
sx q[0];
rz(0.20568307) q[0];
x q[1];
rz(0.63206105) q[2];
sx q[2];
rz(-1.1230317) q[2];
sx q[2];
rz(-2.9438382) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7116075) q[1];
sx q[1];
rz(-1.9031591) q[1];
sx q[1];
rz(-1.3620939) q[1];
rz(-pi) q[2];
rz(-2.9985524) q[3];
sx q[3];
rz(-1.9177799) q[3];
sx q[3];
rz(0.59511371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.70302427) q[2];
sx q[2];
rz(-2.2615304) q[2];
sx q[2];
rz(2.0795889) q[2];
rz(-2.2117173) q[3];
sx q[3];
rz(-2.3558741) q[3];
sx q[3];
rz(0.78833956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36282614) q[0];
sx q[0];
rz(-2.7500948) q[0];
sx q[0];
rz(0.40089259) q[0];
rz(0.76155424) q[1];
sx q[1];
rz(-1.6490033) q[1];
sx q[1];
rz(-2.3525995) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61643314) q[0];
sx q[0];
rz(-3.0130921) q[0];
sx q[0];
rz(-0.12501053) q[0];
rz(-pi) q[1];
rz(1.3775741) q[2];
sx q[2];
rz(-2.4573662) q[2];
sx q[2];
rz(0.88175899) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4915062) q[1];
sx q[1];
rz(-1.6900926) q[1];
sx q[1];
rz(2.5507798) q[1];
rz(-pi) q[2];
rz(-1.5212359) q[3];
sx q[3];
rz(-1.0465099) q[3];
sx q[3];
rz(-1.4567284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.57364982) q[2];
sx q[2];
rz(-1.8486134) q[2];
sx q[2];
rz(0.72081494) q[2];
rz(1.4912262) q[3];
sx q[3];
rz(-2.3590915) q[3];
sx q[3];
rz(1.3097552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1184621) q[0];
sx q[0];
rz(-0.11496249) q[0];
sx q[0];
rz(0.7777099) q[0];
rz(-0.65981162) q[1];
sx q[1];
rz(-0.86645627) q[1];
sx q[1];
rz(-2.7515817) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4068425) q[0];
sx q[0];
rz(-1.5336541) q[0];
sx q[0];
rz(-1.7057306) q[0];
rz(0.74441461) q[2];
sx q[2];
rz(-1.0420024) q[2];
sx q[2];
rz(2.2802558) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.240814) q[1];
sx q[1];
rz(-1.1692746) q[1];
sx q[1];
rz(-2.8177849) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94160752) q[3];
sx q[3];
rz(-2.3857085) q[3];
sx q[3];
rz(-0.98996491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.68221349) q[2];
sx q[2];
rz(-2.7028658) q[2];
sx q[2];
rz(-2.6518346) q[2];
rz(-0.30719906) q[3];
sx q[3];
rz(-0.92370737) q[3];
sx q[3];
rz(-1.3908305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3601892) q[0];
sx q[0];
rz(-1.6506945) q[0];
sx q[0];
rz(-1.6682464) q[0];
rz(2.0669943) q[1];
sx q[1];
rz(-0.85292024) q[1];
sx q[1];
rz(-1.9180752) q[1];
rz(1.9509857) q[2];
sx q[2];
rz(-1.9440704) q[2];
sx q[2];
rz(0.32344641) q[2];
rz(-2.4018039) q[3];
sx q[3];
rz(-0.7480015) q[3];
sx q[3];
rz(1.9960777) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
