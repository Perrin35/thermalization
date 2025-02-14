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
rz(0.30000609) q[0];
sx q[0];
rz(-1.0862779) q[0];
sx q[0];
rz(2.3334184) q[0];
rz(1.9975245) q[1];
sx q[1];
rz(3.9132325) q[1];
sx q[1];
rz(11.003718) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8573541) q[0];
sx q[0];
rz(-0.94737091) q[0];
sx q[0];
rz(2.4592072) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3035266) q[2];
sx q[2];
rz(-1.4216982) q[2];
sx q[2];
rz(-1.9097415) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3951125) q[1];
sx q[1];
rz(-2.0762788) q[1];
sx q[1];
rz(-3.014747) q[1];
rz(-pi) q[2];
rz(-0.64227958) q[3];
sx q[3];
rz(-0.056996973) q[3];
sx q[3];
rz(1.9798724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.943104) q[2];
sx q[2];
rz(-2.9297332) q[2];
sx q[2];
rz(-1.0342106) q[2];
rz(0.69283038) q[3];
sx q[3];
rz(-1.0537909) q[3];
sx q[3];
rz(1.0606162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8387872) q[0];
sx q[0];
rz(-0.028258709) q[0];
sx q[0];
rz(-0.57648188) q[0];
rz(-0.020596404) q[1];
sx q[1];
rz(-2.7437904) q[1];
sx q[1];
rz(1.0284665) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45555025) q[0];
sx q[0];
rz(-0.33838135) q[0];
sx q[0];
rz(0.69940059) q[0];
rz(-pi) q[1];
rz(2.2510301) q[2];
sx q[2];
rz(-0.46698585) q[2];
sx q[2];
rz(2.525223) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0513269) q[1];
sx q[1];
rz(-1.9724066) q[1];
sx q[1];
rz(0.36220596) q[1];
rz(-pi) q[2];
rz(0.50039496) q[3];
sx q[3];
rz(-0.71128363) q[3];
sx q[3];
rz(-0.82373133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.91745201) q[2];
sx q[2];
rz(-2.1805306) q[2];
sx q[2];
rz(1.8769598) q[2];
rz(-0.97359109) q[3];
sx q[3];
rz(-1.6907938) q[3];
sx q[3];
rz(0.26315954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6384386) q[0];
sx q[0];
rz(-1.3042903) q[0];
sx q[0];
rz(0.85357443) q[0];
rz(2.0897934) q[1];
sx q[1];
rz(-2.167326) q[1];
sx q[1];
rz(-0.015017088) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5511709) q[0];
sx q[0];
rz(-2.430924) q[0];
sx q[0];
rz(-1.5842313) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21785801) q[2];
sx q[2];
rz(-3.0431722) q[2];
sx q[2];
rz(-2.276536) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.58850059) q[1];
sx q[1];
rz(-2.4238844) q[1];
sx q[1];
rz(-0.11228447) q[1];
rz(-pi) q[2];
rz(-0.28983966) q[3];
sx q[3];
rz(-1.6320758) q[3];
sx q[3];
rz(2.0743528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0217648) q[2];
sx q[2];
rz(-1.4889577) q[2];
sx q[2];
rz(-2.8601698) q[2];
rz(1.4540539) q[3];
sx q[3];
rz(-1.931793) q[3];
sx q[3];
rz(0.76930261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5948831) q[0];
sx q[0];
rz(-1.2458206) q[0];
sx q[0];
rz(0.17387867) q[0];
rz(-1.3315382) q[1];
sx q[1];
rz(-1.7200108) q[1];
sx q[1];
rz(1.4283659) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4177307) q[0];
sx q[0];
rz(-1.763146) q[0];
sx q[0];
rz(-0.014561166) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7816781) q[2];
sx q[2];
rz(-1.6891306) q[2];
sx q[2];
rz(1.182076) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.62828583) q[1];
sx q[1];
rz(-0.9137972) q[1];
sx q[1];
rz(1.6194067) q[1];
rz(-1.1826199) q[3];
sx q[3];
rz(-1.9132764) q[3];
sx q[3];
rz(1.4701209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4398769) q[2];
sx q[2];
rz(-2.3150847) q[2];
sx q[2];
rz(-0.41713777) q[2];
rz(-2.9719628) q[3];
sx q[3];
rz(-0.093234213) q[3];
sx q[3];
rz(2.4197742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4193264) q[0];
sx q[0];
rz(-0.37501431) q[0];
sx q[0];
rz(-2.8322423) q[0];
rz(0.7695235) q[1];
sx q[1];
rz(-1.0898432) q[1];
sx q[1];
rz(-1.5187029) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2941336) q[0];
sx q[0];
rz(-1.8222229) q[0];
sx q[0];
rz(0.88788106) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3974472) q[2];
sx q[2];
rz(-0.89140201) q[2];
sx q[2];
rz(-1.8216004) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7243821) q[1];
sx q[1];
rz(-1.7599748) q[1];
sx q[1];
rz(-0.030102109) q[1];
rz(-0.937014) q[3];
sx q[3];
rz(-2.975836) q[3];
sx q[3];
rz(2.9637333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7421444) q[2];
sx q[2];
rz(-0.98230201) q[2];
sx q[2];
rz(0.99545109) q[2];
rz(-0.71850145) q[3];
sx q[3];
rz(-1.3281497) q[3];
sx q[3];
rz(2.3464581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28498483) q[0];
sx q[0];
rz(-1.6866848) q[0];
sx q[0];
rz(1.5307776) q[0];
rz(-1.4467622) q[1];
sx q[1];
rz(-1.4434573) q[1];
sx q[1];
rz(1.4379427) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7890678) q[0];
sx q[0];
rz(-1.97856) q[0];
sx q[0];
rz(0.022313281) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9447359) q[2];
sx q[2];
rz(-1.0122006) q[2];
sx q[2];
rz(0.93511673) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.48541486) q[1];
sx q[1];
rz(-0.5159157) q[1];
sx q[1];
rz(-0.1702711) q[1];
rz(-pi) q[2];
rz(-3.1385057) q[3];
sx q[3];
rz(-2.3820138) q[3];
sx q[3];
rz(-2.742033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1630254) q[2];
sx q[2];
rz(-1.7820396) q[2];
sx q[2];
rz(1.2320409) q[2];
rz(1.9949404) q[3];
sx q[3];
rz(-1.8012643) q[3];
sx q[3];
rz(-2.8310217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50800407) q[0];
sx q[0];
rz(-0.81729832) q[0];
sx q[0];
rz(2.001413) q[0];
rz(-2.2445402) q[1];
sx q[1];
rz(-1.3373809) q[1];
sx q[1];
rz(-0.030489771) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59647467) q[0];
sx q[0];
rz(-1.4545355) q[0];
sx q[0];
rz(-1.7932939) q[0];
rz(2.7159436) q[2];
sx q[2];
rz(-2.4535745) q[2];
sx q[2];
rz(-1.4378215) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1709939) q[1];
sx q[1];
rz(-2.355651) q[1];
sx q[1];
rz(1.1517609) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3014841) q[3];
sx q[3];
rz(-1.9193773) q[3];
sx q[3];
rz(-0.19053647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9474779) q[2];
sx q[2];
rz(-2.0674481) q[2];
sx q[2];
rz(-1.2987632) q[2];
rz(-1.4878368) q[3];
sx q[3];
rz(-1.0349118) q[3];
sx q[3];
rz(-1.3207159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79539913) q[0];
sx q[0];
rz(-2.8454056) q[0];
sx q[0];
rz(1.3854223) q[0];
rz(-1.9253383) q[1];
sx q[1];
rz(-1.6938035) q[1];
sx q[1];
rz(-0.48929712) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3922244) q[0];
sx q[0];
rz(-0.67182344) q[0];
sx q[0];
rz(-0.031731204) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1050503) q[2];
sx q[2];
rz(-1.2347691) q[2];
sx q[2];
rz(2.450409) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.41114901) q[1];
sx q[1];
rz(-0.48161067) q[1];
sx q[1];
rz(-1.3661307) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8200666) q[3];
sx q[3];
rz(-1.6695654) q[3];
sx q[3];
rz(2.6840212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3936002) q[2];
sx q[2];
rz(-1.9990498) q[2];
sx q[2];
rz(2.6753814) q[2];
rz(-1.5318058) q[3];
sx q[3];
rz(-1.5038871) q[3];
sx q[3];
rz(-0.32061779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4894678) q[0];
sx q[0];
rz(-2.2418699) q[0];
sx q[0];
rz(-0.049064431) q[0];
rz(-2.9009254) q[1];
sx q[1];
rz(-1.0180232) q[1];
sx q[1];
rz(-3.1191471) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81240679) q[0];
sx q[0];
rz(-2.5617449) q[0];
sx q[0];
rz(-1.6145228) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7540665) q[2];
sx q[2];
rz(-1.2467313) q[2];
sx q[2];
rz(-1.9112196) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6091303) q[1];
sx q[1];
rz(-0.8317308) q[1];
sx q[1];
rz(2.4799281) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4565036) q[3];
sx q[3];
rz(-1.2108608) q[3];
sx q[3];
rz(-1.0459096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0040794) q[2];
sx q[2];
rz(-1.9065607) q[2];
sx q[2];
rz(2.1732886) q[2];
rz(0.83640313) q[3];
sx q[3];
rz(-2.8439549) q[3];
sx q[3];
rz(1.3692726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4150998) q[0];
sx q[0];
rz(-1.2905755) q[0];
sx q[0];
rz(2.7628164) q[0];
rz(3.1355766) q[1];
sx q[1];
rz(-2.6117987) q[1];
sx q[1];
rz(2.5078497) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1155891) q[0];
sx q[0];
rz(-0.77583122) q[0];
sx q[0];
rz(2.6494527) q[0];
rz(-pi) q[1];
x q[1];
rz(1.815755) q[2];
sx q[2];
rz(-2.5542521) q[2];
sx q[2];
rz(0.4191242) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4162916) q[1];
sx q[1];
rz(-1.0868727) q[1];
sx q[1];
rz(2.121622) q[1];
rz(-pi) q[2];
rz(-2.2692436) q[3];
sx q[3];
rz(-2.3309338) q[3];
sx q[3];
rz(-1.8332421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.70574957) q[2];
sx q[2];
rz(-1.5219995) q[2];
sx q[2];
rz(2.5283234) q[2];
rz(-2.4663726) q[3];
sx q[3];
rz(-2.7628511) q[3];
sx q[3];
rz(-0.81775445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2262065) q[0];
sx q[0];
rz(-1.7788667) q[0];
sx q[0];
rz(-1.9052196) q[0];
rz(-2.9551103) q[1];
sx q[1];
rz(-1.3874556) q[1];
sx q[1];
rz(1.5127771) q[1];
rz(1.5507422) q[2];
sx q[2];
rz(-2.0303844) q[2];
sx q[2];
rz(-2.5302237) q[2];
rz(-2.1929838) q[3];
sx q[3];
rz(-1.6453679) q[3];
sx q[3];
rz(-2.9352376) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
