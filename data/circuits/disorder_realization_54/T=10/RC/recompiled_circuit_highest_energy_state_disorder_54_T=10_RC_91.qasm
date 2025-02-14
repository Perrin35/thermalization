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
rz(2.430001) q[0];
sx q[0];
rz(-1.1049668) q[0];
sx q[0];
rz(2.0018863) q[0];
rz(-3.9330339) q[1];
sx q[1];
rz(7.3242261) q[1];
sx q[1];
rz(8.2952226) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72973824) q[0];
sx q[0];
rz(-1.751873) q[0];
sx q[0];
rz(-2.2921261) q[0];
rz(3.1055345) q[2];
sx q[2];
rz(-0.66345045) q[2];
sx q[2];
rz(-2.5306338) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0186386) q[1];
sx q[1];
rz(-1.8644511) q[1];
sx q[1];
rz(-2.3112626) q[1];
rz(-1.6057062) q[3];
sx q[3];
rz(-1.9034317) q[3];
sx q[3];
rz(1.8647609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0396042) q[2];
sx q[2];
rz(-2.3309989) q[2];
sx q[2];
rz(-2.2339036) q[2];
rz(0.11861079) q[3];
sx q[3];
rz(-1.6018931) q[3];
sx q[3];
rz(2.8906726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1949961) q[0];
sx q[0];
rz(-0.66221607) q[0];
sx q[0];
rz(0.10435852) q[0];
rz(-1.9508427) q[1];
sx q[1];
rz(-1.2079116) q[1];
sx q[1];
rz(-0.45713919) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3973244) q[0];
sx q[0];
rz(-1.5400104) q[0];
sx q[0];
rz(-1.7643614) q[0];
x q[1];
rz(-0.89294101) q[2];
sx q[2];
rz(-1.9364076) q[2];
sx q[2];
rz(1.9302238) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3981084) q[1];
sx q[1];
rz(-2.8116075) q[1];
sx q[1];
rz(-1.2742001) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90216402) q[3];
sx q[3];
rz(-1.6794101) q[3];
sx q[3];
rz(1.7028374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.55598688) q[2];
sx q[2];
rz(-1.9827236) q[2];
sx q[2];
rz(-2.9846094) q[2];
rz(0.5213151) q[3];
sx q[3];
rz(-2.3020404) q[3];
sx q[3];
rz(3.1275911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74612015) q[0];
sx q[0];
rz(-2.5465901) q[0];
sx q[0];
rz(0.34573063) q[0];
rz(1.6861457) q[1];
sx q[1];
rz(-1.8691749) q[1];
sx q[1];
rz(1.5706496) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3216305) q[0];
sx q[0];
rz(-1.4819711) q[0];
sx q[0];
rz(1.0246683) q[0];
rz(1.5758031) q[2];
sx q[2];
rz(-0.62267674) q[2];
sx q[2];
rz(1.1345991) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.45528835) q[1];
sx q[1];
rz(-1.521811) q[1];
sx q[1];
rz(0.21360417) q[1];
x q[2];
rz(2.219958) q[3];
sx q[3];
rz(-1.4086257) q[3];
sx q[3];
rz(-2.0258486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7741144) q[2];
sx q[2];
rz(-0.228129) q[2];
sx q[2];
rz(-0.95743123) q[2];
rz(-1.1227603) q[3];
sx q[3];
rz(-1.2343531) q[3];
sx q[3];
rz(-1.7694337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6289309) q[0];
sx q[0];
rz(-1.7191732) q[0];
sx q[0];
rz(-1.5439532) q[0];
rz(-1.7515901) q[1];
sx q[1];
rz(-1.8252239) q[1];
sx q[1];
rz(-1.66473) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8771956) q[0];
sx q[0];
rz(-1.0478643) q[0];
sx q[0];
rz(1.2448798) q[0];
x q[1];
rz(-3.0576287) q[2];
sx q[2];
rz(-1.8083085) q[2];
sx q[2];
rz(-2.824914) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7668666) q[1];
sx q[1];
rz(-1.0488759) q[1];
sx q[1];
rz(0.38066997) q[1];
rz(2.8119726) q[3];
sx q[3];
rz(-1.8330036) q[3];
sx q[3];
rz(-0.58673687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5433898) q[2];
sx q[2];
rz(-1.8566088) q[2];
sx q[2];
rz(-0.99679917) q[2];
rz(1.2654842) q[3];
sx q[3];
rz(-2.4235453) q[3];
sx q[3];
rz(-0.27749458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0604414) q[0];
sx q[0];
rz(-1.5716946) q[0];
sx q[0];
rz(2.0695709) q[0];
rz(-2.187166) q[1];
sx q[1];
rz(-1.9642893) q[1];
sx q[1];
rz(1.6486453) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1562441) q[0];
sx q[0];
rz(-1.3748504) q[0];
sx q[0];
rz(0.77134706) q[0];
x q[1];
rz(-2.6482539) q[2];
sx q[2];
rz(-0.78560053) q[2];
sx q[2];
rz(1.0935022) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6449738) q[1];
sx q[1];
rz(-2.7826834) q[1];
sx q[1];
rz(2.7254057) q[1];
rz(-0.29376438) q[3];
sx q[3];
rz(-2.8076594) q[3];
sx q[3];
rz(1.5042083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91168779) q[2];
sx q[2];
rz(-0.3868843) q[2];
sx q[2];
rz(-1.9913199) q[2];
rz(-3.1005499) q[3];
sx q[3];
rz(-1.5464455) q[3];
sx q[3];
rz(0.40157792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3732442) q[0];
sx q[0];
rz(-2.0337489) q[0];
sx q[0];
rz(-1.9819697) q[0];
rz(1.6805964) q[1];
sx q[1];
rz(-2.9407839) q[1];
sx q[1];
rz(1.2332835) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7338184) q[0];
sx q[0];
rz(-2.3253157) q[0];
sx q[0];
rz(-0.4305779) q[0];
rz(-pi) q[1];
rz(0.057282863) q[2];
sx q[2];
rz(-1.5174688) q[2];
sx q[2];
rz(1.6595226) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6813876) q[1];
sx q[1];
rz(-0.80599313) q[1];
sx q[1];
rz(-2.5780923) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6150675) q[3];
sx q[3];
rz(-1.7987781) q[3];
sx q[3];
rz(2.990852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.636574) q[2];
sx q[2];
rz(-0.36618149) q[2];
sx q[2];
rz(-1.9912857) q[2];
rz(-0.73927528) q[3];
sx q[3];
rz(-1.8940247) q[3];
sx q[3];
rz(0.44481835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66733852) q[0];
sx q[0];
rz(-2.4485454) q[0];
sx q[0];
rz(-0.44878238) q[0];
rz(-1.5397286) q[1];
sx q[1];
rz(-2.0346784) q[1];
sx q[1];
rz(-1.1610228) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13793531) q[0];
sx q[0];
rz(-1.9788392) q[0];
sx q[0];
rz(1.8613226) q[0];
x q[1];
rz(2.6426605) q[2];
sx q[2];
rz(-1.1321486) q[2];
sx q[2];
rz(-2.2427151) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4515069) q[1];
sx q[1];
rz(-1.8645608) q[1];
sx q[1];
rz(2.7965332) q[1];
rz(-pi) q[2];
x q[2];
rz(0.26412873) q[3];
sx q[3];
rz(-2.1089411) q[3];
sx q[3];
rz(1.9028185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1686958) q[2];
sx q[2];
rz(-1.3316414) q[2];
sx q[2];
rz(1.1807582) q[2];
rz(-1.0817889) q[3];
sx q[3];
rz(-1.5321782) q[3];
sx q[3];
rz(-0.42657524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54843724) q[0];
sx q[0];
rz(-1.3189545) q[0];
sx q[0];
rz(0.68761188) q[0];
rz(1.9718735) q[1];
sx q[1];
rz(-0.97427383) q[1];
sx q[1];
rz(2.7489472) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3340281) q[0];
sx q[0];
rz(-1.0962209) q[0];
sx q[0];
rz(-2.9827098) q[0];
rz(-pi) q[1];
rz(0.31198172) q[2];
sx q[2];
rz(-1.7034966) q[2];
sx q[2];
rz(0.080094425) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.029812977) q[1];
sx q[1];
rz(-1.6583484) q[1];
sx q[1];
rz(2.665641) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7358549) q[3];
sx q[3];
rz(-2.9364481) q[3];
sx q[3];
rz(0.9838549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9795867) q[2];
sx q[2];
rz(-2.1350828) q[2];
sx q[2];
rz(1.4957734) q[2];
rz(-1.5227854) q[3];
sx q[3];
rz(-2.3706172) q[3];
sx q[3];
rz(2.5405367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9252121) q[0];
sx q[0];
rz(-0.38327152) q[0];
sx q[0];
rz(-1.3339169) q[0];
rz(1.1434309) q[1];
sx q[1];
rz(-2.3107078) q[1];
sx q[1];
rz(-0.7935895) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.586602) q[0];
sx q[0];
rz(-1.412801) q[0];
sx q[0];
rz(1.642307) q[0];
rz(-1.4758395) q[2];
sx q[2];
rz(-2.5296202) q[2];
sx q[2];
rz(0.31699917) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5836583) q[1];
sx q[1];
rz(-2.0730632) q[1];
sx q[1];
rz(-1.7947547) q[1];
rz(-1.7904441) q[3];
sx q[3];
rz(-2.6815412) q[3];
sx q[3];
rz(-1.3212412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6952343) q[2];
sx q[2];
rz(-2.0909205) q[2];
sx q[2];
rz(-1.0515155) q[2];
rz(2.4255883) q[3];
sx q[3];
rz(-0.99596888) q[3];
sx q[3];
rz(0.21014617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42651549) q[0];
sx q[0];
rz(-2.2861013) q[0];
sx q[0];
rz(0.18095428) q[0];
rz(1.1894233) q[1];
sx q[1];
rz(-1.9633429) q[1];
sx q[1];
rz(-0.98035556) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6342148) q[0];
sx q[0];
rz(-1.6537154) q[0];
sx q[0];
rz(0.50237992) q[0];
rz(-pi) q[1];
x q[1];
rz(0.45578477) q[2];
sx q[2];
rz(-1.96473) q[2];
sx q[2];
rz(1.822871) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7661383) q[1];
sx q[1];
rz(-2.3050637) q[1];
sx q[1];
rz(-1.3260654) q[1];
rz(-0.47360955) q[3];
sx q[3];
rz(-2.7033812) q[3];
sx q[3];
rz(1.1058863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.02562) q[2];
sx q[2];
rz(-0.16416922) q[2];
sx q[2];
rz(1.8170961) q[2];
rz(0.65274158) q[3];
sx q[3];
rz(-2.1721811) q[3];
sx q[3];
rz(0.038221922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41435913) q[0];
sx q[0];
rz(-1.7200732) q[0];
sx q[0];
rz(2.1678069) q[0];
rz(0.7242135) q[1];
sx q[1];
rz(-1.0754633) q[1];
sx q[1];
rz(0.16574688) q[1];
rz(-3.0510256) q[2];
sx q[2];
rz(-1.9048077) q[2];
sx q[2];
rz(2.122369) q[2];
rz(1.6809088) q[3];
sx q[3];
rz(-0.80169818) q[3];
sx q[3];
rz(2.0924951) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
