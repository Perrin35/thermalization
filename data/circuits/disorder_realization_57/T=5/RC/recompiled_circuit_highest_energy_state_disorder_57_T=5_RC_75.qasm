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
rz(-1.3402101) q[0];
sx q[0];
rz(-2.8647381) q[0];
sx q[0];
rz(1.0387596) q[0];
rz(2.7998595) q[1];
sx q[1];
rz(-0.83655292) q[1];
sx q[1];
rz(-0.41681448) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78535372) q[0];
sx q[0];
rz(-0.84945852) q[0];
sx q[0];
rz(-1.4291886) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7567467) q[2];
sx q[2];
rz(-0.59078465) q[2];
sx q[2];
rz(2.7123775) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8326679) q[1];
sx q[1];
rz(-2.1410258) q[1];
sx q[1];
rz(0.69242386) q[1];
rz(1.6031715) q[3];
sx q[3];
rz(-1.7322632) q[3];
sx q[3];
rz(0.77154713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1699528) q[2];
sx q[2];
rz(-1.6728741) q[2];
sx q[2];
rz(1.2539585) q[2];
rz(1.5597255) q[3];
sx q[3];
rz(-0.73837787) q[3];
sx q[3];
rz(-2.019465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9480243) q[0];
sx q[0];
rz(-1.7727611) q[0];
sx q[0];
rz(2.3728306) q[0];
rz(1.3668775) q[1];
sx q[1];
rz(-1.8194852) q[1];
sx q[1];
rz(-2.1034525) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067613451) q[0];
sx q[0];
rz(-1.5609976) q[0];
sx q[0];
rz(-0.0086179535) q[0];
x q[1];
rz(1.9351472) q[2];
sx q[2];
rz(-0.88308217) q[2];
sx q[2];
rz(-1.5503413) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2778216) q[1];
sx q[1];
rz(-1.9367332) q[1];
sx q[1];
rz(0.65524958) q[1];
rz(-pi) q[2];
rz(1.5350545) q[3];
sx q[3];
rz(-0.48887353) q[3];
sx q[3];
rz(1.8349748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.64608964) q[2];
sx q[2];
rz(-1.2868519) q[2];
sx q[2];
rz(0.82026473) q[2];
rz(2.8881554) q[3];
sx q[3];
rz(-2.7866252) q[3];
sx q[3];
rz(-0.36111116) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8609817) q[0];
sx q[0];
rz(-0.51787037) q[0];
sx q[0];
rz(2.2542727) q[0];
rz(-0.53030983) q[1];
sx q[1];
rz(-0.92875004) q[1];
sx q[1];
rz(1.1393772) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6761918) q[0];
sx q[0];
rz(-1.4290819) q[0];
sx q[0];
rz(0.34872524) q[0];
rz(-pi) q[1];
rz(0.33311756) q[2];
sx q[2];
rz(-1.0651759) q[2];
sx q[2];
rz(0.62084711) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.70717057) q[1];
sx q[1];
rz(-1.556273) q[1];
sx q[1];
rz(2.8959031) q[1];
rz(-1.5185131) q[3];
sx q[3];
rz(-1.0526592) q[3];
sx q[3];
rz(0.62925807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.016971074) q[2];
sx q[2];
rz(-1.6861702) q[2];
sx q[2];
rz(-0.43928453) q[2];
rz(-0.47075054) q[3];
sx q[3];
rz(-3.0622523) q[3];
sx q[3];
rz(-1.9973756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73862326) q[0];
sx q[0];
rz(-2.2082177) q[0];
sx q[0];
rz(3.0261107) q[0];
rz(3.0237517) q[1];
sx q[1];
rz(-1.203048) q[1];
sx q[1];
rz(-1.3062564) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8354315) q[0];
sx q[0];
rz(-1.8580313) q[0];
sx q[0];
rz(2.0067418) q[0];
rz(-pi) q[1];
rz(2.4815062) q[2];
sx q[2];
rz(-1.266511) q[2];
sx q[2];
rz(3.0344998) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0769121) q[1];
sx q[1];
rz(-2.2222328) q[1];
sx q[1];
rz(-0.76131911) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4545069) q[3];
sx q[3];
rz(-0.57580417) q[3];
sx q[3];
rz(2.7888014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8370886) q[2];
sx q[2];
rz(-1.3136761) q[2];
sx q[2];
rz(-0.15466776) q[2];
rz(-1.539544) q[3];
sx q[3];
rz(-1.8013026) q[3];
sx q[3];
rz(0.60607564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5681169) q[0];
sx q[0];
rz(-1.0433759) q[0];
sx q[0];
rz(0.83531761) q[0];
rz(2.4256445) q[1];
sx q[1];
rz(-2.7138111) q[1];
sx q[1];
rz(-3.0221525) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9855921) q[0];
sx q[0];
rz(-1.553926) q[0];
sx q[0];
rz(1.5149679) q[0];
x q[1];
rz(0.44942707) q[2];
sx q[2];
rz(-1.7088505) q[2];
sx q[2];
rz(-1.4269478) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2970256) q[1];
sx q[1];
rz(-0.78999619) q[1];
sx q[1];
rz(-0.67761919) q[1];
x q[2];
rz(1.9914845) q[3];
sx q[3];
rz(-1.0771695) q[3];
sx q[3];
rz(2.8378262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9014088) q[2];
sx q[2];
rz(-2.8808424) q[2];
sx q[2];
rz(-0.060997941) q[2];
rz(3.0933464) q[3];
sx q[3];
rz(-1.2333074) q[3];
sx q[3];
rz(2.4129996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7399087) q[0];
sx q[0];
rz(-2.4329199) q[0];
sx q[0];
rz(2.0106864) q[0];
rz(0.4785969) q[1];
sx q[1];
rz(-0.60239783) q[1];
sx q[1];
rz(1.7344249) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76001747) q[0];
sx q[0];
rz(-1.5692466) q[0];
sx q[0];
rz(-0.053561915) q[0];
rz(-pi) q[1];
rz(2.2931104) q[2];
sx q[2];
rz(-1.3568078) q[2];
sx q[2];
rz(-0.011836476) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4987405) q[1];
sx q[1];
rz(-1.7987218) q[1];
sx q[1];
rz(-3.0251316) q[1];
rz(-pi) q[2];
rz(0.70194093) q[3];
sx q[3];
rz(-2.8656368) q[3];
sx q[3];
rz(0.70358932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5799134) q[2];
sx q[2];
rz(-0.40739569) q[2];
sx q[2];
rz(-1.178406) q[2];
rz(1.0493578) q[3];
sx q[3];
rz(-2.6047843) q[3];
sx q[3];
rz(1.8572846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9310164) q[0];
sx q[0];
rz(-1.1897621) q[0];
sx q[0];
rz(-1.806102) q[0];
rz(1.8203075) q[1];
sx q[1];
rz(-1.0711461) q[1];
sx q[1];
rz(2.8335422) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086723344) q[0];
sx q[0];
rz(-0.44741524) q[0];
sx q[0];
rz(2.6347876) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5045583) q[2];
sx q[2];
rz(-0.61372988) q[2];
sx q[2];
rz(1.6047275) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.25702204) q[1];
sx q[1];
rz(-2.0557457) q[1];
sx q[1];
rz(-2.229175) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0940476) q[3];
sx q[3];
rz(-1.1851839) q[3];
sx q[3];
rz(-1.4908718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.41219741) q[2];
sx q[2];
rz(-2.2008379) q[2];
sx q[2];
rz(3.0103053) q[2];
rz(2.4042551) q[3];
sx q[3];
rz(-1.8359343) q[3];
sx q[3];
rz(0.15784119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8312296) q[0];
sx q[0];
rz(-2.0518301) q[0];
sx q[0];
rz(2.6112153) q[0];
rz(1.7165548) q[1];
sx q[1];
rz(-1.1385463) q[1];
sx q[1];
rz(2.4200965) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58396858) q[0];
sx q[0];
rz(-2.1707105) q[0];
sx q[0];
rz(-2.1198089) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4230096) q[2];
sx q[2];
rz(-2.035454) q[2];
sx q[2];
rz(-2.8140659) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4911982) q[1];
sx q[1];
rz(-1.3547239) q[1];
sx q[1];
rz(2.9465527) q[1];
rz(-pi) q[2];
rz(2.3590478) q[3];
sx q[3];
rz(-1.492332) q[3];
sx q[3];
rz(-0.162938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.43712744) q[2];
sx q[2];
rz(-2.7361054) q[2];
sx q[2];
rz(-1.6638157) q[2];
rz(-1.9780698) q[3];
sx q[3];
rz(-2.0587557) q[3];
sx q[3];
rz(1.5014974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8116233) q[0];
sx q[0];
rz(-1.7003308) q[0];
sx q[0];
rz(2.5323618) q[0];
rz(-1.5513264) q[1];
sx q[1];
rz(-0.32569277) q[1];
sx q[1];
rz(-1.4564266) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.698899) q[0];
sx q[0];
rz(-0.44166587) q[0];
sx q[0];
rz(-2.1864258) q[0];
x q[1];
rz(-0.48522075) q[2];
sx q[2];
rz(-1.5169797) q[2];
sx q[2];
rz(-0.14419989) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.63467866) q[1];
sx q[1];
rz(-1.4089481) q[1];
sx q[1];
rz(-1.1062201) q[1];
rz(-2.2942469) q[3];
sx q[3];
rz(-1.0718126) q[3];
sx q[3];
rz(-0.33667281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7234601) q[2];
sx q[2];
rz(-2.2432566) q[2];
sx q[2];
rz(0.88031236) q[2];
rz(-0.20600016) q[3];
sx q[3];
rz(-0.42492953) q[3];
sx q[3];
rz(-2.6515085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3708165) q[0];
sx q[0];
rz(-0.66148615) q[0];
sx q[0];
rz(0.01195512) q[0];
rz(-1.5097584) q[1];
sx q[1];
rz(-0.44523528) q[1];
sx q[1];
rz(-0.59648046) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84688127) q[0];
sx q[0];
rz(-0.46976006) q[0];
sx q[0];
rz(-3.0994979) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5310703) q[2];
sx q[2];
rz(-0.43402616) q[2];
sx q[2];
rz(0.69930017) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47190753) q[1];
sx q[1];
rz(-2.0457256) q[1];
sx q[1];
rz(3.0529939) q[1];
x q[2];
rz(-1.4130728) q[3];
sx q[3];
rz(-1.2765358) q[3];
sx q[3];
rz(-1.0666447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.163588) q[2];
sx q[2];
rz(-2.1808193) q[2];
sx q[2];
rz(2.9445924) q[2];
rz(-1.152285) q[3];
sx q[3];
rz(-0.14324337) q[3];
sx q[3];
rz(-2.6939189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86354179) q[0];
sx q[0];
rz(-2.1320237) q[0];
sx q[0];
rz(-0.19436819) q[0];
rz(-0.20881431) q[1];
sx q[1];
rz(-1.5369692) q[1];
sx q[1];
rz(2.1280638) q[1];
rz(2.2603358) q[2];
sx q[2];
rz(-1.4980346) q[2];
sx q[2];
rz(-3.0537506) q[2];
rz(-0.98255929) q[3];
sx q[3];
rz(-2.0053902) q[3];
sx q[3];
rz(-2.7872661) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
