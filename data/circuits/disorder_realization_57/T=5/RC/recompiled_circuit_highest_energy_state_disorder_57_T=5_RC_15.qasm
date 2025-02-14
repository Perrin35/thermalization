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
rz(2.7998595) q[1];
sx q[1];
rz(-0.83655292) q[1];
sx q[1];
rz(-0.41681448) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57274103) q[0];
sx q[0];
rz(-2.4089455) q[0];
sx q[0];
rz(-0.1591263) q[0];
rz(-pi) q[1];
rz(-2.5854163) q[2];
sx q[2];
rz(-1.7814629) q[2];
sx q[2];
rz(-0.81708252) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3089247) q[1];
sx q[1];
rz(-2.1410258) q[1];
sx q[1];
rz(0.69242386) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9454222) q[3];
sx q[3];
rz(-2.9769398) q[3];
sx q[3];
rz(0.97033721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9716399) q[2];
sx q[2];
rz(-1.4687186) q[2];
sx q[2];
rz(-1.2539585) q[2];
rz(-1.5818671) q[3];
sx q[3];
rz(-0.73837787) q[3];
sx q[3];
rz(-2.019465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1935683) q[0];
sx q[0];
rz(-1.7727611) q[0];
sx q[0];
rz(2.3728306) q[0];
rz(1.7747152) q[1];
sx q[1];
rz(-1.3221075) q[1];
sx q[1];
rz(1.0381402) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067613451) q[0];
sx q[0];
rz(-1.5609976) q[0];
sx q[0];
rz(-0.0086179535) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2064455) q[2];
sx q[2];
rz(-0.88308217) q[2];
sx q[2];
rz(-1.5912513) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27127344) q[1];
sx q[1];
rz(-2.4044577) q[1];
sx q[1];
rz(-0.56136521) q[1];
rz(-pi) q[2];
rz(1.0821876) q[3];
sx q[3];
rz(-1.587579) q[3];
sx q[3];
rz(-0.29573655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.495503) q[2];
sx q[2];
rz(-1.2868519) q[2];
sx q[2];
rz(0.82026473) q[2];
rz(0.25343728) q[3];
sx q[3];
rz(-2.7866252) q[3];
sx q[3];
rz(0.36111116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28061098) q[0];
sx q[0];
rz(-2.6237223) q[0];
sx q[0];
rz(0.88731998) q[0];
rz(-2.6112828) q[1];
sx q[1];
rz(-2.2128426) q[1];
sx q[1];
rz(-2.0022154) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6761918) q[0];
sx q[0];
rz(-1.4290819) q[0];
sx q[0];
rz(-0.34872524) q[0];
x q[1];
rz(-2.8084751) q[2];
sx q[2];
rz(-1.0651759) q[2];
sx q[2];
rz(0.62084711) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2816086) q[1];
sx q[1];
rz(-1.8164595) q[1];
sx q[1];
rz(1.5857693) q[1];
rz(1.5185131) q[3];
sx q[3];
rz(-2.0889335) q[3];
sx q[3];
rz(0.62925807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1246216) q[2];
sx q[2];
rz(-1.4554224) q[2];
sx q[2];
rz(-2.7023081) q[2];
rz(-0.47075054) q[3];
sx q[3];
rz(-0.079340383) q[3];
sx q[3];
rz(-1.144217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4029694) q[0];
sx q[0];
rz(-0.93337494) q[0];
sx q[0];
rz(-3.0261107) q[0];
rz(-0.11784095) q[1];
sx q[1];
rz(-1.9385447) q[1];
sx q[1];
rz(-1.8353362) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4233526) q[0];
sx q[0];
rz(-0.51694345) q[0];
sx q[0];
rz(0.96036185) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9491862) q[2];
sx q[2];
rz(-2.1956964) q[2];
sx q[2];
rz(-1.4493799) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0688248) q[1];
sx q[1];
rz(-0.95736527) q[1];
sx q[1];
rz(-2.3062506) q[1];
rz(-pi) q[2];
rz(-0.99808295) q[3];
sx q[3];
rz(-1.6340165) q[3];
sx q[3];
rz(-1.8259189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8370886) q[2];
sx q[2];
rz(-1.8279165) q[2];
sx q[2];
rz(-0.15466776) q[2];
rz(-1.539544) q[3];
sx q[3];
rz(-1.3402901) q[3];
sx q[3];
rz(-0.60607564) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5681169) q[0];
sx q[0];
rz(-1.0433759) q[0];
sx q[0];
rz(-2.306275) q[0];
rz(0.71594816) q[1];
sx q[1];
rz(-0.42778152) q[1];
sx q[1];
rz(-3.0221525) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41573856) q[0];
sx q[0];
rz(-1.5149759) q[0];
sx q[0];
rz(0.016896642) q[0];
x q[1];
rz(-2.6921656) q[2];
sx q[2];
rz(-1.7088505) q[2];
sx q[2];
rz(-1.4269478) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9000098) q[1];
sx q[1];
rz(-1.1092343) q[1];
sx q[1];
rz(2.475283) q[1];
rz(1.1501081) q[3];
sx q[3];
rz(-1.0771695) q[3];
sx q[3];
rz(0.30376645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9014088) q[2];
sx q[2];
rz(-2.8808424) q[2];
sx q[2];
rz(3.0805947) q[2];
rz(0.048246233) q[3];
sx q[3];
rz(-1.2333074) q[3];
sx q[3];
rz(0.72859305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7399087) q[0];
sx q[0];
rz(-0.70867276) q[0];
sx q[0];
rz(-2.0106864) q[0];
rz(0.4785969) q[1];
sx q[1];
rz(-0.60239783) q[1];
sx q[1];
rz(1.7344249) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3307307) q[0];
sx q[0];
rz(-1.5172345) q[0];
sx q[0];
rz(-1.5723482) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8883775) q[2];
sx q[2];
rz(-0.7478315) q[2];
sx q[2];
rz(-1.795447) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9764912) q[1];
sx q[1];
rz(-0.25549421) q[1];
sx q[1];
rz(2.0352023) q[1];
rz(-pi) q[2];
rz(-0.21295548) q[3];
sx q[3];
rz(-1.3939438) q[3];
sx q[3];
rz(2.9573754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.56167928) q[2];
sx q[2];
rz(-2.734197) q[2];
sx q[2];
rz(1.9631867) q[2];
rz(-2.0922349) q[3];
sx q[3];
rz(-0.5368084) q[3];
sx q[3];
rz(1.2843081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2105763) q[0];
sx q[0];
rz(-1.9518305) q[0];
sx q[0];
rz(-1.3354906) q[0];
rz(1.3212851) q[1];
sx q[1];
rz(-1.0711461) q[1];
sx q[1];
rz(-2.8335422) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0548693) q[0];
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
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8561473) q[1];
sx q[1];
rz(-0.79570192) q[1];
sx q[1];
rz(-2.281762) q[1];
rz(1.9568029) q[3];
sx q[3];
rz(-1.6148477) q[3];
sx q[3];
rz(0.06202997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7293952) q[2];
sx q[2];
rz(-0.9407548) q[2];
sx q[2];
rz(3.0103053) q[2];
rz(-0.73733759) q[3];
sx q[3];
rz(-1.8359343) q[3];
sx q[3];
rz(-2.9837515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8312296) q[0];
sx q[0];
rz(-2.0518301) q[0];
sx q[0];
rz(-0.53037733) q[0];
rz(-1.4250379) q[1];
sx q[1];
rz(-1.1385463) q[1];
sx q[1];
rz(-0.72149611) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3193766) q[0];
sx q[0];
rz(-1.1255029) q[0];
sx q[0];
rz(-0.67586835) q[0];
rz(2.6725298) q[2];
sx q[2];
rz(-1.7028168) q[2];
sx q[2];
rz(1.309883) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.65039448) q[1];
sx q[1];
rz(-1.3547239) q[1];
sx q[1];
rz(2.9465527) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0305392) q[3];
sx q[3];
rz(-0.78563443) q[3];
sx q[3];
rz(1.3291886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.43712744) q[2];
sx q[2];
rz(-0.4054873) q[2];
sx q[2];
rz(-1.4777769) q[2];
rz(-1.9780698) q[3];
sx q[3];
rz(-2.0587557) q[3];
sx q[3];
rz(1.5014974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32996938) q[0];
sx q[0];
rz(-1.7003308) q[0];
sx q[0];
rz(2.5323618) q[0];
rz(-1.5902663) q[1];
sx q[1];
rz(-2.8158999) q[1];
sx q[1];
rz(-1.4564266) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1066172) q[0];
sx q[0];
rz(-1.2143232) q[0];
sx q[0];
rz(-2.8750505) q[0];
rz(2.6563719) q[2];
sx q[2];
rz(-1.5169797) q[2];
sx q[2];
rz(-0.14419989) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2472154) q[1];
sx q[1];
rz(-0.4900107) q[1];
sx q[1];
rz(1.2213329) q[1];
rz(-pi) q[2];
rz(0.62868613) q[3];
sx q[3];
rz(-0.95052023) q[3];
sx q[3];
rz(1.6339906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.41813254) q[2];
sx q[2];
rz(-0.89833608) q[2];
sx q[2];
rz(2.2612803) q[2];
rz(-0.20600016) q[3];
sx q[3];
rz(-0.42492953) q[3];
sx q[3];
rz(0.4900842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.77077615) q[0];
sx q[0];
rz(-0.66148615) q[0];
sx q[0];
rz(3.1296375) q[0];
rz(1.6318343) q[1];
sx q[1];
rz(-2.6963574) q[1];
sx q[1];
rz(0.59648046) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79967989) q[0];
sx q[0];
rz(-2.0401067) q[0];
sx q[0];
rz(-1.5494359) q[0];
rz(1.311074) q[2];
sx q[2];
rz(-1.2190281) q[2];
sx q[2];
rz(0.042366926) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.66377992) q[1];
sx q[1];
rz(-2.6590909) q[1];
sx q[1];
rz(-1.4003808) q[1];
rz(0.29774547) q[3];
sx q[3];
rz(-1.4199054) q[3];
sx q[3];
rz(2.6835364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.163588) q[2];
sx q[2];
rz(-2.1808193) q[2];
sx q[2];
rz(2.9445924) q[2];
rz(1.152285) q[3];
sx q[3];
rz(-0.14324337) q[3];
sx q[3];
rz(2.6939189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86354179) q[0];
sx q[0];
rz(-2.1320237) q[0];
sx q[0];
rz(-0.19436819) q[0];
rz(2.9327783) q[1];
sx q[1];
rz(-1.5369692) q[1];
sx q[1];
rz(2.1280638) q[1];
rz(0.094194407) q[2];
sx q[2];
rz(-2.258156) q[2];
sx q[2];
rz(-1.5428262) q[2];
rz(-2.267425) q[3];
sx q[3];
rz(-0.71577358) q[3];
sx q[3];
rz(-1.7795455) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
