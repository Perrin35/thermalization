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
rz(1.6765321) q[0];
sx q[0];
rz(-2.9002011) q[0];
sx q[0];
rz(-0.13225947) q[0];
rz(-0.013068696) q[1];
sx q[1];
rz(-0.71813923) q[1];
sx q[1];
rz(3.1225966) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4723139) q[0];
sx q[0];
rz(-1.3558071) q[0];
sx q[0];
rz(-0.90256067) q[0];
rz(-pi) q[1];
rz(-1.1055345) q[2];
sx q[2];
rz(-2.9334515) q[2];
sx q[2];
rz(-0.15424745) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.63848313) q[1];
sx q[1];
rz(-1.3967525) q[1];
sx q[1];
rz(2.8139504) q[1];
rz(-pi) q[2];
rz(-0.98145841) q[3];
sx q[3];
rz(-0.27772003) q[3];
sx q[3];
rz(-0.026175682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0211109) q[2];
sx q[2];
rz(-2.3349473) q[2];
sx q[2];
rz(2.3497904) q[2];
rz(1.0857438) q[3];
sx q[3];
rz(-1.9224242) q[3];
sx q[3];
rz(-2.048548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2708112) q[0];
sx q[0];
rz(-2.9830611) q[0];
sx q[0];
rz(-0.32919163) q[0];
rz(1.5022494) q[1];
sx q[1];
rz(-2.2396125) q[1];
sx q[1];
rz(0.42218581) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49445186) q[0];
sx q[0];
rz(-1.1472771) q[0];
sx q[0];
rz(0.5943055) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5936826) q[2];
sx q[2];
rz(-1.413125) q[2];
sx q[2];
rz(-1.6239177) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.364155) q[1];
sx q[1];
rz(-1.6171439) q[1];
sx q[1];
rz(-2.8175161) q[1];
rz(-0.77236891) q[3];
sx q[3];
rz(-1.6447146) q[3];
sx q[3];
rz(-2.071601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4210356) q[2];
sx q[2];
rz(-1.1853508) q[2];
sx q[2];
rz(1.1400918) q[2];
rz(0.0044048443) q[3];
sx q[3];
rz(-1.5811698) q[3];
sx q[3];
rz(-3.0018023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36188257) q[0];
sx q[0];
rz(-3.0411868) q[0];
sx q[0];
rz(2.7959339) q[0];
rz(-2.0480305) q[1];
sx q[1];
rz(-0.82229096) q[1];
sx q[1];
rz(2.0872769) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7525268) q[0];
sx q[0];
rz(-1.9574165) q[0];
sx q[0];
rz(1.2091314) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81297154) q[2];
sx q[2];
rz(-1.6457498) q[2];
sx q[2];
rz(-3.093442) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7400764) q[1];
sx q[1];
rz(-1.1187727) q[1];
sx q[1];
rz(-1.6416618) q[1];
rz(-0.90462138) q[3];
sx q[3];
rz(-1.5890749) q[3];
sx q[3];
rz(-1.8005023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.67636079) q[2];
sx q[2];
rz(-0.54764843) q[2];
sx q[2];
rz(-2.5993627) q[2];
rz(-0.17255653) q[3];
sx q[3];
rz(-1.6514643) q[3];
sx q[3];
rz(-2.0358613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3365823) q[0];
sx q[0];
rz(-1.7886826) q[0];
sx q[0];
rz(-1.8804869) q[0];
rz(2.005596) q[1];
sx q[1];
rz(-2.0295862) q[1];
sx q[1];
rz(-0.13066185) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.703637) q[0];
sx q[0];
rz(-1.2017631) q[0];
sx q[0];
rz(-2.8519408) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50640743) q[2];
sx q[2];
rz(-1.7676438) q[2];
sx q[2];
rz(-1.0723237) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.16176341) q[1];
sx q[1];
rz(-1.5391304) q[1];
sx q[1];
rz(2.0791441) q[1];
rz(-pi) q[2];
rz(2.4774136) q[3];
sx q[3];
rz(-1.2746433) q[3];
sx q[3];
rz(0.68063762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.38264349) q[2];
sx q[2];
rz(-3.068277) q[2];
sx q[2];
rz(0.22155133) q[2];
rz(2.4394636) q[3];
sx q[3];
rz(-1.5626855) q[3];
sx q[3];
rz(-2.4041972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5704983) q[0];
sx q[0];
rz(-1.8966738) q[0];
sx q[0];
rz(0.13394295) q[0];
rz(0.36930034) q[1];
sx q[1];
rz(-1.1993273) q[1];
sx q[1];
rz(-1.7514924) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1655859) q[0];
sx q[0];
rz(-2.0780141) q[0];
sx q[0];
rz(2.8602703) q[0];
rz(-2.8763387) q[2];
sx q[2];
rz(-1.2758288) q[2];
sx q[2];
rz(-0.8910999) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.304313) q[1];
sx q[1];
rz(-1.2943448) q[1];
sx q[1];
rz(2.5778092) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69306121) q[3];
sx q[3];
rz(-0.45940889) q[3];
sx q[3];
rz(1.0074248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.2964581) q[2];
sx q[2];
rz(-0.90201169) q[2];
sx q[2];
rz(-1.585539) q[2];
rz(-2.8499917) q[3];
sx q[3];
rz(-0.4762989) q[3];
sx q[3];
rz(-0.27879032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45048243) q[0];
sx q[0];
rz(-0.99523681) q[0];
sx q[0];
rz(0.52571785) q[0];
rz(-2.763343) q[1];
sx q[1];
rz(-1.5317761) q[1];
sx q[1];
rz(-2.8674616) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7129242) q[0];
sx q[0];
rz(-0.98085058) q[0];
sx q[0];
rz(1.5398067) q[0];
x q[1];
rz(2.7080665) q[2];
sx q[2];
rz(-1.9872403) q[2];
sx q[2];
rz(0.64420494) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3124481) q[1];
sx q[1];
rz(-0.48407468) q[1];
sx q[1];
rz(0.99906724) q[1];
x q[2];
rz(1.0438221) q[3];
sx q[3];
rz(-0.29114215) q[3];
sx q[3];
rz(-3.1028735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5452925) q[2];
sx q[2];
rz(-1.2240852) q[2];
sx q[2];
rz(-2.8833585) q[2];
rz(-1.5457414) q[3];
sx q[3];
rz(-1.3629379) q[3];
sx q[3];
rz(0.405092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3426568) q[0];
sx q[0];
rz(-0.5181784) q[0];
sx q[0];
rz(-0.10547353) q[0];
rz(2.8666829) q[1];
sx q[1];
rz(-1.7173488) q[1];
sx q[1];
rz(-2.8555433) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8409922) q[0];
sx q[0];
rz(-1.4541885) q[0];
sx q[0];
rz(-0.18516115) q[0];
rz(-2.3069165) q[2];
sx q[2];
rz(-1.3501985) q[2];
sx q[2];
rz(-0.59360628) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.090324245) q[1];
sx q[1];
rz(-0.39127884) q[1];
sx q[1];
rz(1.4011739) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5887542) q[3];
sx q[3];
rz(-1.1854608) q[3];
sx q[3];
rz(2.5401153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.99453551) q[2];
sx q[2];
rz(-2.4592082) q[2];
sx q[2];
rz(2.6915754) q[2];
rz(-0.58722812) q[3];
sx q[3];
rz(-1.8000032) q[3];
sx q[3];
rz(-1.2500866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29379544) q[0];
sx q[0];
rz(-1.6814517) q[0];
sx q[0];
rz(1.9548804) q[0];
rz(-2.8222491) q[1];
sx q[1];
rz(-1.0217383) q[1];
sx q[1];
rz(0.26225463) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38459331) q[0];
sx q[0];
rz(-2.5525064) q[0];
sx q[0];
rz(-0.37269816) q[0];
rz(-pi) q[1];
rz(-1.9840711) q[2];
sx q[2];
rz(-1.6954633) q[2];
sx q[2];
rz(1.1050129) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.1848534) q[1];
sx q[1];
rz(-1.0295086) q[1];
sx q[1];
rz(-0.095603099) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3877034) q[3];
sx q[3];
rz(-2.5544832) q[3];
sx q[3];
rz(1.0894437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.35683262) q[2];
sx q[2];
rz(-2.6748952) q[2];
sx q[2];
rz(-0.32877767) q[2];
rz(-1.6262936) q[3];
sx q[3];
rz(-1.3582151) q[3];
sx q[3];
rz(2.7330107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3560155) q[0];
sx q[0];
rz(-2.6308036) q[0];
sx q[0];
rz(-0.27895862) q[0];
rz(-0.51697671) q[1];
sx q[1];
rz(-0.37716436) q[1];
sx q[1];
rz(1.7163716) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12869975) q[0];
sx q[0];
rz(-0.94966799) q[0];
sx q[0];
rz(1.6245882) q[0];
rz(0.89377706) q[2];
sx q[2];
rz(-0.89134673) q[2];
sx q[2];
rz(0.51365863) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.60630262) q[1];
sx q[1];
rz(-0.53724506) q[1];
sx q[1];
rz(0.37588889) q[1];
rz(-pi) q[2];
rz(0.29028671) q[3];
sx q[3];
rz(-2.2072574) q[3];
sx q[3];
rz(0.94062128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7384537) q[2];
sx q[2];
rz(-2.8454915) q[2];
sx q[2];
rz(-0.1599172) q[2];
rz(-2.3219409) q[3];
sx q[3];
rz(-1.9194226) q[3];
sx q[3];
rz(-0.73614365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4495471) q[0];
sx q[0];
rz(-1.8501546) q[0];
sx q[0];
rz(0.29577574) q[0];
rz(2.6462818) q[1];
sx q[1];
rz(-0.43423978) q[1];
sx q[1];
rz(-1.6326509) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9072278) q[0];
sx q[0];
rz(-1.9669269) q[0];
sx q[0];
rz(-2.7559126) q[0];
rz(-pi) q[1];
rz(1.2184754) q[2];
sx q[2];
rz(-1.5218456) q[2];
sx q[2];
rz(1.5201598) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.26824046) q[1];
sx q[1];
rz(-0.11246364) q[1];
sx q[1];
rz(0.64551533) q[1];
x q[2];
rz(0.4949244) q[3];
sx q[3];
rz(-1.9175745) q[3];
sx q[3];
rz(-1.0241672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8136924) q[2];
sx q[2];
rz(-1.4699961) q[2];
sx q[2];
rz(0.63391614) q[2];
rz(3.0564803) q[3];
sx q[3];
rz(-1.8952993) q[3];
sx q[3];
rz(2.3688721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2278628) q[0];
sx q[0];
rz(-1.5416523) q[0];
sx q[0];
rz(-0.1097485) q[0];
rz(1.9204503) q[1];
sx q[1];
rz(-0.84538645) q[1];
sx q[1];
rz(0.56199817) q[1];
rz(-2.0094677) q[2];
sx q[2];
rz(-1.3668458) q[2];
sx q[2];
rz(0.20076164) q[2];
rz(3.1172995) q[3];
sx q[3];
rz(-1.3391936) q[3];
sx q[3];
rz(-1.5486123) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
