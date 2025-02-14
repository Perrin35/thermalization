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
rz(-2.3251301) q[0];
sx q[0];
rz(-0.10179585) q[0];
sx q[0];
rz(2.6020004) q[0];
rz(0.49630961) q[1];
sx q[1];
rz(-0.30975431) q[1];
sx q[1];
rz(0.53909477) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64083034) q[0];
sx q[0];
rz(-2.0120562) q[0];
sx q[0];
rz(1.5426969) q[0];
x q[1];
rz(-1.4194054) q[2];
sx q[2];
rz(-1.2149723) q[2];
sx q[2];
rz(2.9265938) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7254759) q[1];
sx q[1];
rz(-1.2486287) q[1];
sx q[1];
rz(2.9941818) q[1];
x q[2];
rz(-3.059637) q[3];
sx q[3];
rz(-2.2181619) q[3];
sx q[3];
rz(-2.820197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.76217905) q[2];
sx q[2];
rz(-0.87612408) q[2];
sx q[2];
rz(-1.1890821) q[2];
rz(-1.1842229) q[3];
sx q[3];
rz(-0.90457478) q[3];
sx q[3];
rz(1.5466461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9963843) q[0];
sx q[0];
rz(-1.5518016) q[0];
sx q[0];
rz(-1.8183964) q[0];
rz(-2.6595751) q[1];
sx q[1];
rz(-0.91164416) q[1];
sx q[1];
rz(-0.97420305) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8708987) q[0];
sx q[0];
rz(-0.93241954) q[0];
sx q[0];
rz(0.22298546) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0423099) q[2];
sx q[2];
rz(-1.6834604) q[2];
sx q[2];
rz(0.4140062) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.99813733) q[1];
sx q[1];
rz(-2.1051198) q[1];
sx q[1];
rz(1.2972531) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.810435) q[3];
sx q[3];
rz(-1.0175034) q[3];
sx q[3];
rz(0.026913337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.0049858967) q[2];
sx q[2];
rz(-0.58711457) q[2];
sx q[2];
rz(0.74964398) q[2];
rz(2.7096115) q[3];
sx q[3];
rz(-2.0260729) q[3];
sx q[3];
rz(1.8384793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5169446) q[0];
sx q[0];
rz(-0.2722781) q[0];
sx q[0];
rz(-2.5352449) q[0];
rz(1.3336522) q[1];
sx q[1];
rz(-1.9275815) q[1];
sx q[1];
rz(-2.8820754) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1071651) q[0];
sx q[0];
rz(-1.3239064) q[0];
sx q[0];
rz(-2.9811923) q[0];
x q[1];
rz(0.86028966) q[2];
sx q[2];
rz(-0.26488556) q[2];
sx q[2];
rz(0.67240326) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1705679) q[1];
sx q[1];
rz(-1.1491421) q[1];
sx q[1];
rz(0.41280156) q[1];
rz(-pi) q[2];
rz(-0.93575017) q[3];
sx q[3];
rz(-0.66422909) q[3];
sx q[3];
rz(-2.2838044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2567265) q[2];
sx q[2];
rz(-1.3639516) q[2];
sx q[2];
rz(1.5941031) q[2];
rz(0.78440845) q[3];
sx q[3];
rz(-1.514785) q[3];
sx q[3];
rz(-1.5659531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6467658) q[0];
sx q[0];
rz(-2.0443199) q[0];
sx q[0];
rz(-0.25949091) q[0];
rz(1.2207813) q[1];
sx q[1];
rz(-2.6916598) q[1];
sx q[1];
rz(1.4422013) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83373678) q[0];
sx q[0];
rz(-2.4015732) q[0];
sx q[0];
rz(0.83777512) q[0];
x q[1];
rz(-3.0109826) q[2];
sx q[2];
rz(-0.97960661) q[2];
sx q[2];
rz(1.3863877) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71953785) q[1];
sx q[1];
rz(-1.2201078) q[1];
sx q[1];
rz(1.3077875) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2529834) q[3];
sx q[3];
rz(-1.5164135) q[3];
sx q[3];
rz(2.3171484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0356902) q[2];
sx q[2];
rz(-1.8041958) q[2];
sx q[2];
rz(-2.7316459) q[2];
rz(-1.9299054) q[3];
sx q[3];
rz(-2.4544921) q[3];
sx q[3];
rz(-1.3338026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.425151) q[0];
sx q[0];
rz(-2.4254159) q[0];
sx q[0];
rz(0.27405611) q[0];
rz(0.43824276) q[1];
sx q[1];
rz(-1.2934877) q[1];
sx q[1];
rz(-0.73572198) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5925795) q[0];
sx q[0];
rz(-1.4837449) q[0];
sx q[0];
rz(-0.75496952) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76916285) q[2];
sx q[2];
rz(-0.20649466) q[2];
sx q[2];
rz(2.1392876) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2949038) q[1];
sx q[1];
rz(-1.3219993) q[1];
sx q[1];
rz(1.2685246) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51208074) q[3];
sx q[3];
rz(-2.1554073) q[3];
sx q[3];
rz(1.1325815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2030486) q[2];
sx q[2];
rz(-1.4578578) q[2];
sx q[2];
rz(-2.5277444) q[2];
rz(-2.7326873) q[3];
sx q[3];
rz(-0.96858612) q[3];
sx q[3];
rz(0.74234211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.3530389) q[0];
sx q[0];
rz(-0.63816324) q[0];
sx q[0];
rz(1.9045389) q[0];
rz(-1.6963814) q[1];
sx q[1];
rz(-2.684869) q[1];
sx q[1];
rz(2.9663185) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5492903) q[0];
sx q[0];
rz(-1.3978551) q[0];
sx q[0];
rz(-0.025385234) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66368033) q[2];
sx q[2];
rz(-2.366407) q[2];
sx q[2];
rz(0.41942393) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3172463) q[1];
sx q[1];
rz(-0.69686705) q[1];
sx q[1];
rz(2.4319355) q[1];
x q[2];
rz(-0.57717429) q[3];
sx q[3];
rz(-1.0374616) q[3];
sx q[3];
rz(1.2791025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1767629) q[2];
sx q[2];
rz(-1.333933) q[2];
sx q[2];
rz(-2.1691587) q[2];
rz(-1.6644299) q[3];
sx q[3];
rz(-1.9921314) q[3];
sx q[3];
rz(0.76178637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2870188) q[0];
sx q[0];
rz(-0.022495689) q[0];
sx q[0];
rz(-0.35312411) q[0];
rz(1.0278206) q[1];
sx q[1];
rz(-1.9408344) q[1];
sx q[1];
rz(-1.0110528) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7142657) q[0];
sx q[0];
rz(-1.2135226) q[0];
sx q[0];
rz(1.5861804) q[0];
rz(-0.89247668) q[2];
sx q[2];
rz(-0.8304285) q[2];
sx q[2];
rz(2.3811946) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.87349975) q[1];
sx q[1];
rz(-0.38485369) q[1];
sx q[1];
rz(1.9995688) q[1];
x q[2];
rz(1.2561463) q[3];
sx q[3];
rz(-1.8487747) q[3];
sx q[3];
rz(-1.6772456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.31356835) q[2];
sx q[2];
rz(-0.32005969) q[2];
sx q[2];
rz(1.774452) q[2];
rz(1.8732871) q[3];
sx q[3];
rz(-1.4788078) q[3];
sx q[3];
rz(-0.84793276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1139514) q[0];
sx q[0];
rz(-2.5340762) q[0];
sx q[0];
rz(1.129958) q[0];
rz(0.21895151) q[1];
sx q[1];
rz(-1.7349225) q[1];
sx q[1];
rz(-2.2311282) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2856871) q[0];
sx q[0];
rz(-0.98074161) q[0];
sx q[0];
rz(1.9737712) q[0];
rz(-2.9772308) q[2];
sx q[2];
rz(-2.6577762) q[2];
sx q[2];
rz(-1.9828895) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7607019) q[1];
sx q[1];
rz(-2.5538951) q[1];
sx q[1];
rz(-2.0668849) q[1];
rz(-pi) q[2];
rz(-2.4164356) q[3];
sx q[3];
rz(-1.2598383) q[3];
sx q[3];
rz(2.6127315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3160481) q[2];
sx q[2];
rz(-0.80222183) q[2];
sx q[2];
rz(2.3160589) q[2];
rz(-0.46142203) q[3];
sx q[3];
rz(-1.2666707) q[3];
sx q[3];
rz(2.6958444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5017186) q[0];
sx q[0];
rz(-1.8032782) q[0];
sx q[0];
rz(0.83183944) q[0];
rz(-2.4141451) q[1];
sx q[1];
rz(-1.7214382) q[1];
sx q[1];
rz(-0.096253455) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98672685) q[0];
sx q[0];
rz(-2.4270202) q[0];
sx q[0];
rz(2.3113219) q[0];
rz(-2.9092838) q[2];
sx q[2];
rz(-2.8806928) q[2];
sx q[2];
rz(2.1292854) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.61273328) q[1];
sx q[1];
rz(-1.3599571) q[1];
sx q[1];
rz(-0.092540578) q[1];
rz(-1.7670125) q[3];
sx q[3];
rz(-1.726578) q[3];
sx q[3];
rz(2.2957612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7353797) q[2];
sx q[2];
rz(-1.7071743) q[2];
sx q[2];
rz(-1.5489138) q[2];
rz(2.5088572) q[3];
sx q[3];
rz(-1.9097208) q[3];
sx q[3];
rz(-0.24937853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7525472) q[0];
sx q[0];
rz(-1.0405552) q[0];
sx q[0];
rz(2.7885875) q[0];
rz(-2.8189335) q[1];
sx q[1];
rz(-0.32538515) q[1];
sx q[1];
rz(-0.37193146) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4034525) q[0];
sx q[0];
rz(-2.7785025) q[0];
sx q[0];
rz(-0.9842666) q[0];
rz(-1.7499128) q[2];
sx q[2];
rz(-1.2460016) q[2];
sx q[2];
rz(-1.3472404) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0363716) q[1];
sx q[1];
rz(-0.37215713) q[1];
sx q[1];
rz(2.6879265) q[1];
rz(-pi) q[2];
rz(0.84623611) q[3];
sx q[3];
rz(-1.5414943) q[3];
sx q[3];
rz(-1.8054363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9670664) q[2];
sx q[2];
rz(-1.1897949) q[2];
sx q[2];
rz(1.8444427) q[2];
rz(2.2938812) q[3];
sx q[3];
rz(-1.984237) q[3];
sx q[3];
rz(2.1367836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90947718) q[0];
sx q[0];
rz(-1.7625325) q[0];
sx q[0];
rz(0.43176227) q[0];
rz(-1.7535946) q[1];
sx q[1];
rz(-1.9276062) q[1];
sx q[1];
rz(3.066317) q[1];
rz(2.4405932) q[2];
sx q[2];
rz(-0.95078118) q[2];
sx q[2];
rz(-2.4581428) q[2];
rz(0.41107117) q[3];
sx q[3];
rz(-2.3007994) q[3];
sx q[3];
rz(2.3074335) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
