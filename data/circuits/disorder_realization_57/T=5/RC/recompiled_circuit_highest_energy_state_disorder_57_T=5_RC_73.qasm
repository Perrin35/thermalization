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
rz(1.8013826) q[0];
sx q[0];
rz(-0.27685452) q[0];
sx q[0];
rz(2.1028331) q[0];
rz(-0.34173319) q[1];
sx q[1];
rz(3.9781456) q[1];
sx q[1];
rz(9.8415924) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2622803) q[0];
sx q[0];
rz(-1.4646155) q[0];
sx q[0];
rz(2.4152629) q[0];
x q[1];
rz(1.3241346) q[2];
sx q[2];
rz(-2.1132872) q[2];
sx q[2];
rz(0.88298029) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.15910251) q[1];
sx q[1];
rz(-1.0034518) q[1];
sx q[1];
rz(-0.87615396) q[1];
rz(-pi) q[2];
rz(-0.19617041) q[3];
sx q[3];
rz(-0.16465287) q[3];
sx q[3];
rz(2.1712554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9716399) q[2];
sx q[2];
rz(-1.4687186) q[2];
sx q[2];
rz(-1.2539585) q[2];
rz(-1.5597255) q[3];
sx q[3];
rz(-2.4032148) q[3];
sx q[3];
rz(1.1221277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1935683) q[0];
sx q[0];
rz(-1.3688315) q[0];
sx q[0];
rz(2.3728306) q[0];
rz(1.7747152) q[1];
sx q[1];
rz(-1.8194852) q[1];
sx q[1];
rz(2.1034525) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5032673) q[0];
sx q[0];
rz(-1.5621788) q[0];
sx q[0];
rz(1.5609972) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7323167) q[2];
sx q[2];
rz(-2.3774494) q[2];
sx q[2];
rz(1.0502358) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.97538725) q[1];
sx q[1];
rz(-0.96549368) q[1];
sx q[1];
rz(1.1206085) q[1];
rz(-pi) q[2];
rz(1.5350545) q[3];
sx q[3];
rz(-0.48887353) q[3];
sx q[3];
rz(-1.3066178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.64608964) q[2];
sx q[2];
rz(-1.2868519) q[2];
sx q[2];
rz(0.82026473) q[2];
rz(-0.25343728) q[3];
sx q[3];
rz(-2.7866252) q[3];
sx q[3];
rz(-0.36111116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28061098) q[0];
sx q[0];
rz(-2.6237223) q[0];
sx q[0];
rz(-0.88731998) q[0];
rz(-0.53030983) q[1];
sx q[1];
rz(-0.92875004) q[1];
sx q[1];
rz(-2.0022154) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.087505) q[0];
sx q[0];
rz(-1.9158792) q[0];
sx q[0];
rz(-1.4201384) q[0];
rz(-2.1007295) q[2];
sx q[2];
rz(-1.2806674) q[2];
sx q[2];
rz(-1.1159971) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.92147747) q[1];
sx q[1];
rz(-0.24610983) q[1];
sx q[1];
rz(-0.059644894) q[1];
x q[2];
rz(0.091413012) q[3];
sx q[3];
rz(-0.52052906) q[3];
sx q[3];
rz(-0.5239858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1246216) q[2];
sx q[2];
rz(-1.4554224) q[2];
sx q[2];
rz(0.43928453) q[2];
rz(-2.6708421) q[3];
sx q[3];
rz(-3.0622523) q[3];
sx q[3];
rz(1.9973756) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73862326) q[0];
sx q[0];
rz(-2.2082177) q[0];
sx q[0];
rz(-3.0261107) q[0];
rz(-3.0237517) q[1];
sx q[1];
rz(-1.9385447) q[1];
sx q[1];
rz(-1.3062564) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30616111) q[0];
sx q[0];
rz(-1.2835614) q[0];
sx q[0];
rz(2.0067418) q[0];
x q[1];
rz(-0.66008644) q[2];
sx q[2];
rz(-1.8750817) q[2];
sx q[2];
rz(-3.0344998) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1237643) q[1];
sx q[1];
rz(-2.1514847) q[1];
sx q[1];
rz(2.3821217) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4545069) q[3];
sx q[3];
rz(-2.5657885) q[3];
sx q[3];
rz(0.35279122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.304504) q[2];
sx q[2];
rz(-1.3136761) q[2];
sx q[2];
rz(2.9869249) q[2];
rz(-1.6020487) q[3];
sx q[3];
rz(-1.8013026) q[3];
sx q[3];
rz(2.535517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5681169) q[0];
sx q[0];
rz(-2.0982168) q[0];
sx q[0];
rz(-0.83531761) q[0];
rz(-2.4256445) q[1];
sx q[1];
rz(-0.42778152) q[1];
sx q[1];
rz(-3.0221525) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41573856) q[0];
sx q[0];
rz(-1.5149759) q[0];
sx q[0];
rz(-0.016896642) q[0];
rz(-pi) q[1];
rz(0.44942707) q[2];
sx q[2];
rz(-1.7088505) q[2];
sx q[2];
rz(1.7146448) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1491974) q[1];
sx q[1];
rz(-0.98434224) q[1];
sx q[1];
rz(-2.1349364) q[1];
rz(0.53262696) q[3];
sx q[3];
rz(-1.9386734) q[3];
sx q[3];
rz(-2.0834578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2401838) q[2];
sx q[2];
rz(-0.26075026) q[2];
sx q[2];
rz(3.0805947) q[2];
rz(3.0933464) q[3];
sx q[3];
rz(-1.9082853) q[3];
sx q[3];
rz(-2.4129996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.401684) q[0];
sx q[0];
rz(-2.4329199) q[0];
sx q[0];
rz(-1.1309062) q[0];
rz(0.4785969) q[1];
sx q[1];
rz(-0.60239783) q[1];
sx q[1];
rz(1.7344249) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78188183) q[0];
sx q[0];
rz(-0.053584307) q[0];
sx q[0];
rz(3.1126541) q[0];
rz(0.84848225) q[2];
sx q[2];
rz(-1.7847848) q[2];
sx q[2];
rz(3.1297562) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.045627198) q[1];
sx q[1];
rz(-1.4573604) q[1];
sx q[1];
rz(1.8002224) q[1];
rz(1.3899441) q[3];
sx q[3];
rz(-1.3612124) q[3];
sx q[3];
rz(1.7169894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5799134) q[2];
sx q[2];
rz(-0.40739569) q[2];
sx q[2];
rz(-1.9631867) q[2];
rz(2.0922349) q[3];
sx q[3];
rz(-0.5368084) q[3];
sx q[3];
rz(-1.2843081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2105763) q[0];
sx q[0];
rz(-1.1897621) q[0];
sx q[0];
rz(1.806102) q[0];
rz(-1.8203075) q[1];
sx q[1];
rz(-2.0704465) q[1];
sx q[1];
rz(2.8335422) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.193452) q[0];
sx q[0];
rz(-1.7823671) q[0];
sx q[0];
rz(-0.39724644) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5045583) q[2];
sx q[2];
rz(-2.5278628) q[2];
sx q[2];
rz(-1.6047275) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.96770331) q[1];
sx q[1];
rz(-2.1427665) q[1];
sx q[1];
rz(2.553945) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0940476) q[3];
sx q[3];
rz(-1.1851839) q[3];
sx q[3];
rz(-1.6507208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7293952) q[2];
sx q[2];
rz(-0.9407548) q[2];
sx q[2];
rz(-3.0103053) q[2];
rz(0.73733759) q[3];
sx q[3];
rz(-1.8359343) q[3];
sx q[3];
rz(-0.15784119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-0.8312296) q[0];
sx q[0];
rz(-2.0518301) q[0];
sx q[0];
rz(-2.6112153) q[0];
rz(1.4250379) q[1];
sx q[1];
rz(-1.1385463) q[1];
sx q[1];
rz(-2.4200965) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24151267) q[0];
sx q[0];
rz(-2.3519313) q[0];
sx q[0];
rz(0.65171839) q[0];
rz(-pi) q[1];
rz(1.4230096) q[2];
sx q[2];
rz(-2.035454) q[2];
sx q[2];
rz(-2.8140659) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7466654) q[1];
sx q[1];
rz(-2.8515366) q[1];
sx q[1];
rz(0.8474838) q[1];
rz(-pi) q[2];
rz(1.6812229) q[3];
sx q[3];
rz(-0.79130606) q[3];
sx q[3];
rz(1.4856389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.43712744) q[2];
sx q[2];
rz(-2.7361054) q[2];
sx q[2];
rz(1.4777769) q[2];
rz(-1.9780698) q[3];
sx q[3];
rz(-2.0587557) q[3];
sx q[3];
rz(1.5014974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8116233) q[0];
sx q[0];
rz(-1.7003308) q[0];
sx q[0];
rz(0.60923088) q[0];
rz(-1.5902663) q[1];
sx q[1];
rz(-0.32569277) q[1];
sx q[1];
rz(-1.685166) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0349755) q[0];
sx q[0];
rz(-1.9272695) q[0];
sx q[0];
rz(2.8750505) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5099735) q[2];
sx q[2];
rz(-2.0552539) q[2];
sx q[2];
rz(1.4549507) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2472154) q[1];
sx q[1];
rz(-0.4900107) q[1];
sx q[1];
rz(1.2213329) q[1];
rz(-2.5129065) q[3];
sx q[3];
rz(-2.1910724) q[3];
sx q[3];
rz(-1.6339906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.41813254) q[2];
sx q[2];
rz(-2.2432566) q[2];
sx q[2];
rz(2.2612803) q[2];
rz(-2.9355925) q[3];
sx q[3];
rz(-0.42492953) q[3];
sx q[3];
rz(-0.4900842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77077615) q[0];
sx q[0];
rz(-2.4801065) q[0];
sx q[0];
rz(-3.1296375) q[0];
rz(1.6318343) q[1];
sx q[1];
rz(-0.44523528) q[1];
sx q[1];
rz(2.5451122) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84688127) q[0];
sx q[0];
rz(-0.46976006) q[0];
sx q[0];
rz(0.042094783) q[0];
rz(-pi) q[1];
rz(-1.311074) q[2];
sx q[2];
rz(-1.2190281) q[2];
sx q[2];
rz(3.0992257) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.47190753) q[1];
sx q[1];
rz(-1.0958671) q[1];
sx q[1];
rz(0.088598786) q[1];
x q[2];
rz(1.4130728) q[3];
sx q[3];
rz(-1.8650569) q[3];
sx q[3];
rz(2.074948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.163588) q[2];
sx q[2];
rz(-0.96077335) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2780509) q[0];
sx q[0];
rz(-2.1320237) q[0];
sx q[0];
rz(-0.19436819) q[0];
rz(-2.9327783) q[1];
sx q[1];
rz(-1.6046235) q[1];
sx q[1];
rz(-1.0135289) q[1];
rz(-1.456719) q[2];
sx q[2];
rz(-0.69274215) q[2];
sx q[2];
rz(-1.3950166) q[2];
rz(0.50894751) q[3];
sx q[3];
rz(-1.043368) q[3];
sx q[3];
rz(2.1989078) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
