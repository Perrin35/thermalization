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
rz(0.99474466) q[0];
sx q[0];
rz(3.7288546) q[0];
sx q[0];
rz(10.050339) q[0];
rz(0.62043959) q[1];
sx q[1];
rz(4.1320463) q[1];
sx q[1];
rz(10.208701) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44255689) q[0];
sx q[0];
rz(-0.43414206) q[0];
sx q[0];
rz(0.0084369466) q[0];
rz(-pi) q[1];
x q[1];
rz(3.002499) q[2];
sx q[2];
rz(-1.4882097) q[2];
sx q[2];
rz(0.54143426) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9770324) q[1];
sx q[1];
rz(-0.98233084) q[1];
sx q[1];
rz(2.2684437) q[1];
x q[2];
rz(-0.77772452) q[3];
sx q[3];
rz(-1.5081662) q[3];
sx q[3];
rz(2.3304813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1938532) q[2];
sx q[2];
rz(-2.096602) q[2];
sx q[2];
rz(2.4742773) q[2];
rz(2.8269178) q[3];
sx q[3];
rz(-0.59942013) q[3];
sx q[3];
rz(2.2144894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84745234) q[0];
sx q[0];
rz(-0.84646928) q[0];
sx q[0];
rz(-0.2997998) q[0];
rz(-1.0046129) q[1];
sx q[1];
rz(-0.71669465) q[1];
sx q[1];
rz(-2.2915548) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4097828) q[0];
sx q[0];
rz(-2.1730521) q[0];
sx q[0];
rz(-0.11645634) q[0];
x q[1];
rz(-1.4910788) q[2];
sx q[2];
rz(-1.4007916) q[2];
sx q[2];
rz(-1.0159279) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0850573) q[1];
sx q[1];
rz(-1.3621482) q[1];
sx q[1];
rz(-1.2056808) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1095457) q[3];
sx q[3];
rz(-0.40478313) q[3];
sx q[3];
rz(-1.7657775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5689508) q[2];
sx q[2];
rz(-2.6093542) q[2];
sx q[2];
rz(2.4786095) q[2];
rz(-0.7524544) q[3];
sx q[3];
rz(-1.3198676) q[3];
sx q[3];
rz(0.19645709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3982584) q[0];
sx q[0];
rz(-0.4011811) q[0];
sx q[0];
rz(2.4617526) q[0];
rz(-0.72194779) q[1];
sx q[1];
rz(-2.5847021) q[1];
sx q[1];
rz(0.78071761) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3762622) q[0];
sx q[0];
rz(-1.1701382) q[0];
sx q[0];
rz(3.1332934) q[0];
rz(-pi) q[1];
rz(3.137792) q[2];
sx q[2];
rz(-1.6569306) q[2];
sx q[2];
rz(-1.116556) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5131849) q[1];
sx q[1];
rz(-0.76090136) q[1];
sx q[1];
rz(-2.582798) q[1];
rz(-pi) q[2];
rz(2.9227681) q[3];
sx q[3];
rz(-2.2940562) q[3];
sx q[3];
rz(-1.7909602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9146933) q[2];
sx q[2];
rz(-1.9136027) q[2];
sx q[2];
rz(0.66960382) q[2];
rz(0.90034825) q[3];
sx q[3];
rz(-1.0122296) q[3];
sx q[3];
rz(-2.3646234) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25903073) q[0];
sx q[0];
rz(-0.43864033) q[0];
sx q[0];
rz(2.2341109) q[0];
rz(-2.5471845) q[1];
sx q[1];
rz(-0.92955697) q[1];
sx q[1];
rz(2.0118735) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8639899) q[0];
sx q[0];
rz(-2.4765795) q[0];
sx q[0];
rz(1.8825085) q[0];
rz(-pi) q[1];
rz(-1.3200892) q[2];
sx q[2];
rz(-0.54566979) q[2];
sx q[2];
rz(0.16599338) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2767375) q[1];
sx q[1];
rz(-1.3001137) q[1];
sx q[1];
rz(2.5321743) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80813415) q[3];
sx q[3];
rz(-2.7459347) q[3];
sx q[3];
rz(1.7803223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.10924673) q[2];
sx q[2];
rz(-1.3999908) q[2];
sx q[2];
rz(0.67231154) q[2];
rz(-0.0191056) q[3];
sx q[3];
rz(-2.8002383) q[3];
sx q[3];
rz(0.13489558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30447176) q[0];
sx q[0];
rz(-0.79455513) q[0];
sx q[0];
rz(2.9735612) q[0];
rz(1.2270323) q[1];
sx q[1];
rz(-1.9214168) q[1];
sx q[1];
rz(-2.8438445) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7452324) q[0];
sx q[0];
rz(-2.0848635) q[0];
sx q[0];
rz(-2.5395509) q[0];
rz(-2.1476168) q[2];
sx q[2];
rz(-1.0712396) q[2];
sx q[2];
rz(0.039488878) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7682957) q[1];
sx q[1];
rz(-1.3650286) q[1];
sx q[1];
rz(1.5962257) q[1];
x q[2];
rz(-1.7295804) q[3];
sx q[3];
rz(-1.5701236) q[3];
sx q[3];
rz(2.1070045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9510522) q[2];
sx q[2];
rz(-1.0995883) q[2];
sx q[2];
rz(0.61868787) q[2];
rz(0.80858532) q[3];
sx q[3];
rz(-2.6638668) q[3];
sx q[3];
rz(-1.3396858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.488778) q[0];
sx q[0];
rz(-1.5050911) q[0];
sx q[0];
rz(-0.46522796) q[0];
rz(-2.6564927) q[1];
sx q[1];
rz(-0.41050375) q[1];
sx q[1];
rz(-0.088265158) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7951435) q[0];
sx q[0];
rz(-0.69899054) q[0];
sx q[0];
rz(-2.3153852) q[0];
rz(-pi) q[1];
rz(1.5563929) q[2];
sx q[2];
rz(-2.4126137) q[2];
sx q[2];
rz(-2.5165503) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3151206) q[1];
sx q[1];
rz(-1.406028) q[1];
sx q[1];
rz(-0.75854782) q[1];
rz(-pi) q[2];
rz(1.4824105) q[3];
sx q[3];
rz(-1.8575577) q[3];
sx q[3];
rz(0.6500611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7469067) q[2];
sx q[2];
rz(-2.2862819) q[2];
sx q[2];
rz(-2.6439903) q[2];
rz(-0.48007128) q[3];
sx q[3];
rz(-2.2205455) q[3];
sx q[3];
rz(-3.0729496) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29480544) q[0];
sx q[0];
rz(-2.6755896) q[0];
sx q[0];
rz(-1.4303327) q[0];
rz(1.826674) q[1];
sx q[1];
rz(-2.7480835) q[1];
sx q[1];
rz(1.5865145) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2908132) q[0];
sx q[0];
rz(-1.4516136) q[0];
sx q[0];
rz(-1.9656154) q[0];
rz(0.58760037) q[2];
sx q[2];
rz(-2.1955904) q[2];
sx q[2];
rz(-1.8163565) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.58007121) q[1];
sx q[1];
rz(-0.59021705) q[1];
sx q[1];
rz(1.3886222) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6756546) q[3];
sx q[3];
rz(-0.63408454) q[3];
sx q[3];
rz(0.88632562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.84307182) q[2];
sx q[2];
rz(-2.4303747) q[2];
sx q[2];
rz(-2.2141875) q[2];
rz(-1.1585506) q[3];
sx q[3];
rz(-2.4535593) q[3];
sx q[3];
rz(-0.096175171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4226828) q[0];
sx q[0];
rz(-0.89253187) q[0];
sx q[0];
rz(2.6655777) q[0];
rz(-2.1503275) q[1];
sx q[1];
rz(-0.74949336) q[1];
sx q[1];
rz(0.083866619) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96362417) q[0];
sx q[0];
rz(-1.2126847) q[0];
sx q[0];
rz(1.0492508) q[0];
rz(1.2071361) q[2];
sx q[2];
rz(-2.0061099) q[2];
sx q[2];
rz(2.2029623) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.36271154) q[1];
sx q[1];
rz(-1.485752) q[1];
sx q[1];
rz(-1.9910732) q[1];
rz(1.7937035) q[3];
sx q[3];
rz(-0.37316868) q[3];
sx q[3];
rz(0.56115967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91741651) q[2];
sx q[2];
rz(-0.26599628) q[2];
sx q[2];
rz(0.64845294) q[2];
rz(0.90483856) q[3];
sx q[3];
rz(-1.4584533) q[3];
sx q[3];
rz(-2.8308433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34173486) q[0];
sx q[0];
rz(-0.71001995) q[0];
sx q[0];
rz(-2.572686) q[0];
rz(-0.83373249) q[1];
sx q[1];
rz(-1.2833475) q[1];
sx q[1];
rz(1.0047097) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89812714) q[0];
sx q[0];
rz(-1.3782189) q[0];
sx q[0];
rz(-2.5054727) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9453853) q[2];
sx q[2];
rz(-0.56972144) q[2];
sx q[2];
rz(-1.2562871) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6697996) q[1];
sx q[1];
rz(-0.87369117) q[1];
sx q[1];
rz(3.0473188) q[1];
x q[2];
rz(-1.6528204) q[3];
sx q[3];
rz(-2.1567705) q[3];
sx q[3];
rz(2.2524633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4363165) q[2];
sx q[2];
rz(-2.269561) q[2];
sx q[2];
rz(-2.7936068) q[2];
rz(-2.3157388) q[3];
sx q[3];
rz(-2.9233942) q[3];
sx q[3];
rz(0.60514778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0843049) q[0];
sx q[0];
rz(-2.075752) q[0];
sx q[0];
rz(-0.2070981) q[0];
rz(-2.6026978) q[1];
sx q[1];
rz(-0.62108827) q[1];
sx q[1];
rz(-0.60212392) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83610632) q[0];
sx q[0];
rz(-2.0417111) q[0];
sx q[0];
rz(-2.6013646) q[0];
x q[1];
rz(2.5670577) q[2];
sx q[2];
rz(-0.67107397) q[2];
sx q[2];
rz(-0.35071638) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.94127266) q[1];
sx q[1];
rz(-1.4956022) q[1];
sx q[1];
rz(-1.978051) q[1];
rz(-pi) q[2];
rz(-0.35134985) q[3];
sx q[3];
rz(-2.3627639) q[3];
sx q[3];
rz(-3.1064649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3567317) q[2];
sx q[2];
rz(-2.1897903) q[2];
sx q[2];
rz(1.3447364) q[2];
rz(2.6209659) q[3];
sx q[3];
rz(-2.8713363) q[3];
sx q[3];
rz(-2.3394913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6945334) q[0];
sx q[0];
rz(-2.4166528) q[0];
sx q[0];
rz(1.7550533) q[0];
rz(-0.78320349) q[1];
sx q[1];
rz(-1.7192817) q[1];
sx q[1];
rz(-1.5789938) q[1];
rz(2.798525) q[2];
sx q[2];
rz(-1.82556) q[2];
sx q[2];
rz(-1.0221046) q[2];
rz(2.9014723) q[3];
sx q[3];
rz(-1.5594395) q[3];
sx q[3];
rz(-1.6394564) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
