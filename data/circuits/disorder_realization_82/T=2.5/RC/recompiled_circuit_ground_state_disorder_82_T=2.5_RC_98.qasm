OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.3759484) q[0];
sx q[0];
rz(-2.7274237) q[0];
sx q[0];
rz(0.8015269) q[0];
rz(3.1385359) q[1];
sx q[1];
rz(-2.2771775) q[1];
sx q[1];
rz(-3.0473696) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.191358) q[0];
sx q[0];
rz(-1.9468029) q[0];
sx q[0];
rz(-2.3467499) q[0];
rz(-pi) q[1];
rz(1.4613737) q[2];
sx q[2];
rz(-2.8673807) q[2];
sx q[2];
rz(2.4502129) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7677418) q[1];
sx q[1];
rz(-0.5597502) q[1];
sx q[1];
rz(1.00882) q[1];
x q[2];
rz(-1.5685969) q[3];
sx q[3];
rz(-1.4181976) q[3];
sx q[3];
rz(-2.8650177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.5376544) q[2];
sx q[2];
rz(-2.5014169) q[2];
sx q[2];
rz(1.630265) q[2];
rz(0.18621914) q[3];
sx q[3];
rz(-0.43034601) q[3];
sx q[3];
rz(-0.65096861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6179825) q[0];
sx q[0];
rz(-1.1815434) q[0];
sx q[0];
rz(-3.004916) q[0];
rz(0.094820529) q[1];
sx q[1];
rz(-2.6780728) q[1];
sx q[1];
rz(-0.76900855) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0083197) q[0];
sx q[0];
rz(-1.1964487) q[0];
sx q[0];
rz(1.4248614) q[0];
rz(-pi) q[1];
rz(1.2142173) q[2];
sx q[2];
rz(-2.2327023) q[2];
sx q[2];
rz(-1.8476768) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.63751672) q[1];
sx q[1];
rz(-2.3055162) q[1];
sx q[1];
rz(1.7903663) q[1];
rz(-1.7058825) q[3];
sx q[3];
rz(-1.5363494) q[3];
sx q[3];
rz(-2.5273539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.39917055) q[2];
sx q[2];
rz(-0.48226446) q[2];
sx q[2];
rz(1.7102309) q[2];
rz(-2.6509905) q[3];
sx q[3];
rz(-0.49121818) q[3];
sx q[3];
rz(-0.77559364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8785777) q[0];
sx q[0];
rz(-0.37913015) q[0];
sx q[0];
rz(2.0468792) q[0];
rz(1.3021944) q[1];
sx q[1];
rz(-0.98919386) q[1];
sx q[1];
rz(-3.1375695) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72633171) q[0];
sx q[0];
rz(-1.7021966) q[0];
sx q[0];
rz(-0.8432998) q[0];
rz(-pi) q[1];
rz(2.268774) q[2];
sx q[2];
rz(-2.1790163) q[2];
sx q[2];
rz(1.5814511) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.074306369) q[1];
sx q[1];
rz(-1.3137215) q[1];
sx q[1];
rz(-1.6021614) q[1];
x q[2];
rz(-2.4278042) q[3];
sx q[3];
rz(-1.7587489) q[3];
sx q[3];
rz(-2.1655708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5502988) q[2];
sx q[2];
rz(-0.58612263) q[2];
sx q[2];
rz(2.6934521) q[2];
rz(-0.48629931) q[3];
sx q[3];
rz(-2.2794162) q[3];
sx q[3];
rz(2.015131) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22911856) q[0];
sx q[0];
rz(-1.4294701) q[0];
sx q[0];
rz(-0.66977704) q[0];
rz(-1.9122596) q[1];
sx q[1];
rz(-0.45893097) q[1];
sx q[1];
rz(2.5965447) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4507283) q[0];
sx q[0];
rz(-1.04869) q[0];
sx q[0];
rz(1.0342811) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82370211) q[2];
sx q[2];
rz(-2.3266226) q[2];
sx q[2];
rz(0.65836775) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5814171) q[1];
sx q[1];
rz(-2.9184074) q[1];
sx q[1];
rz(1.2944004) q[1];
rz(-0.7983851) q[3];
sx q[3];
rz(-1.2252062) q[3];
sx q[3];
rz(-1.8963199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5547319) q[2];
sx q[2];
rz(-1.6394337) q[2];
sx q[2];
rz(2.9736605) q[2];
rz(1.933291) q[3];
sx q[3];
rz(-0.30253634) q[3];
sx q[3];
rz(1.4471853) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5912882) q[0];
sx q[0];
rz(-1.2606324) q[0];
sx q[0];
rz(0.0038797832) q[0];
rz(2.8344391) q[1];
sx q[1];
rz(-2.3852564) q[1];
sx q[1];
rz(-0.53877962) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8545495) q[0];
sx q[0];
rz(-1.7998621) q[0];
sx q[0];
rz(-1.4255217) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8454119) q[2];
sx q[2];
rz(-2.6950172) q[2];
sx q[2];
rz(2.2889529) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5153577) q[1];
sx q[1];
rz(-1.9222676) q[1];
sx q[1];
rz(0.61351794) q[1];
rz(-pi) q[2];
x q[2];
rz(2.410048) q[3];
sx q[3];
rz(-1.9172102) q[3];
sx q[3];
rz(0.15958318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6946081) q[2];
sx q[2];
rz(-2.0422523) q[2];
sx q[2];
rz(1.851932) q[2];
rz(-1.6192216) q[3];
sx q[3];
rz(-0.74519849) q[3];
sx q[3];
rz(-1.4815909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30037844) q[0];
sx q[0];
rz(-2.2024246) q[0];
sx q[0];
rz(-3.0389431) q[0];
rz(-2.4094021) q[1];
sx q[1];
rz(-2.3468572) q[1];
sx q[1];
rz(-0.57714677) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10023663) q[0];
sx q[0];
rz(-1.5994342) q[0];
sx q[0];
rz(-0.0036544637) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4056751) q[2];
sx q[2];
rz(-0.99011865) q[2];
sx q[2];
rz(-0.48027793) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.093955479) q[1];
sx q[1];
rz(-0.13493061) q[1];
sx q[1];
rz(0.62420242) q[1];
rz(0.34586819) q[3];
sx q[3];
rz(-0.78639275) q[3];
sx q[3];
rz(-2.6707471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.010043667) q[2];
sx q[2];
rz(-2.4250344) q[2];
sx q[2];
rz(-1.5232167) q[2];
rz(0.067954436) q[3];
sx q[3];
rz(-2.8165635) q[3];
sx q[3];
rz(0.84097356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0129358) q[0];
sx q[0];
rz(-1.034863) q[0];
sx q[0];
rz(-0.77142429) q[0];
rz(-2.6029288) q[1];
sx q[1];
rz(-1.3464876) q[1];
sx q[1];
rz(1.0661941) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2832269) q[0];
sx q[0];
rz(-1.8233426) q[0];
sx q[0];
rz(0.079312328) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0577002) q[2];
sx q[2];
rz(-1.3670245) q[2];
sx q[2];
rz(0.17901006) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.33246189) q[1];
sx q[1];
rz(-1.5839424) q[1];
sx q[1];
rz(2.9244208) q[1];
x q[2];
rz(0.33447075) q[3];
sx q[3];
rz(-1.3190509) q[3];
sx q[3];
rz(1.9289964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6923339) q[2];
sx q[2];
rz(-2.3350495) q[2];
sx q[2];
rz(2.6991357) q[2];
rz(-2.9509406) q[3];
sx q[3];
rz(-0.72398829) q[3];
sx q[3];
rz(2.382522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8964748) q[0];
sx q[0];
rz(-2.9560095) q[0];
sx q[0];
rz(-1.8316487) q[0];
rz(2.7438121) q[1];
sx q[1];
rz(-0.82292992) q[1];
sx q[1];
rz(1.7479053) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1289541) q[0];
sx q[0];
rz(-1.1348179) q[0];
sx q[0];
rz(-0.97431824) q[0];
rz(-1.2692503) q[2];
sx q[2];
rz(-1.796495) q[2];
sx q[2];
rz(1.7762453) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.85227784) q[1];
sx q[1];
rz(-0.56486579) q[1];
sx q[1];
rz(-0.58128618) q[1];
rz(-pi) q[2];
x q[2];
rz(1.769479) q[3];
sx q[3];
rz(-1.4666838) q[3];
sx q[3];
rz(-0.35752359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.028367793) q[2];
sx q[2];
rz(-0.71030474) q[2];
sx q[2];
rz(-3.0681211) q[2];
rz(-0.69463378) q[3];
sx q[3];
rz(-1.2725384) q[3];
sx q[3];
rz(-2.1139483) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46591127) q[0];
sx q[0];
rz(-0.2008734) q[0];
sx q[0];
rz(0.95941108) q[0];
rz(-2.361946) q[1];
sx q[1];
rz(-2.1387073) q[1];
sx q[1];
rz(-2.8693105) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1355541) q[0];
sx q[0];
rz(-1.5840496) q[0];
sx q[0];
rz(-0.0051965836) q[0];
x q[1];
rz(1.9302888) q[2];
sx q[2];
rz(-0.73795107) q[2];
sx q[2];
rz(2.8518845) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4705369) q[1];
sx q[1];
rz(-1.1644286) q[1];
sx q[1];
rz(2.0996835) q[1];
rz(-pi) q[2];
rz(-0.97314463) q[3];
sx q[3];
rz(-0.91734195) q[3];
sx q[3];
rz(-2.0367095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.73725629) q[2];
sx q[2];
rz(-1.5055089) q[2];
sx q[2];
rz(-0.20757248) q[2];
rz(3.0370144) q[3];
sx q[3];
rz(-2.6360377) q[3];
sx q[3];
rz(1.526621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.680147) q[0];
sx q[0];
rz(-1.2274281) q[0];
sx q[0];
rz(2.0848059) q[0];
rz(-0.54569221) q[1];
sx q[1];
rz(-2.1634407) q[1];
sx q[1];
rz(0.32870764) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0775111) q[0];
sx q[0];
rz(-1.9010549) q[0];
sx q[0];
rz(-1.0868549) q[0];
x q[1];
rz(0.07241197) q[2];
sx q[2];
rz(-1.5161523) q[2];
sx q[2];
rz(-1.3018198) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4663965) q[1];
sx q[1];
rz(-1.5874784) q[1];
sx q[1];
rz(1.2287324) q[1];
x q[2];
rz(-1.8449109) q[3];
sx q[3];
rz(-0.81820993) q[3];
sx q[3];
rz(-0.75358281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9823965) q[2];
sx q[2];
rz(-0.13267645) q[2];
sx q[2];
rz(2.0170434) q[2];
rz(3.0647965) q[3];
sx q[3];
rz(-0.43282893) q[3];
sx q[3];
rz(-0.92397773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15976739) q[0];
sx q[0];
rz(-1.2564909) q[0];
sx q[0];
rz(1.5782574) q[0];
rz(-0.81623296) q[1];
sx q[1];
rz(-1.7385794) q[1];
sx q[1];
rz(2.0508918) q[1];
rz(-0.74538415) q[2];
sx q[2];
rz(-1.1911662) q[2];
sx q[2];
rz(0.62698812) q[2];
rz(-1.5745926) q[3];
sx q[3];
rz(-1.1654614) q[3];
sx q[3];
rz(-0.75335791) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
