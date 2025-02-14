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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1358937) q[0];
sx q[0];
rz(-1.5672475) q[0];
sx q[0];
rz(-2.7074642) q[0];
rz(-2.6033635) q[2];
sx q[2];
rz(-2.9799649) q[2];
sx q[2];
rz(-1.5618351) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9770324) q[1];
sx q[1];
rz(-0.98233084) q[1];
sx q[1];
rz(-2.2684437) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4830071) q[3];
sx q[3];
rz(-0.79500073) q[3];
sx q[3];
rz(-2.3203497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1938532) q[2];
sx q[2];
rz(-1.0449907) q[2];
sx q[2];
rz(-2.4742773) q[2];
rz(2.8269178) q[3];
sx q[3];
rz(-2.5421725) q[3];
sx q[3];
rz(-2.2144894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2941403) q[0];
sx q[0];
rz(-2.2951234) q[0];
sx q[0];
rz(0.2997998) q[0];
rz(-2.1369797) q[1];
sx q[1];
rz(-0.71669465) q[1];
sx q[1];
rz(2.2915548) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4097828) q[0];
sx q[0];
rz(-2.1730521) q[0];
sx q[0];
rz(3.0251363) q[0];
rz(2.9710567) q[2];
sx q[2];
rz(-1.6493622) q[2];
sx q[2];
rz(2.6002392) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5932754) q[1];
sx q[1];
rz(-1.9276345) q[1];
sx q[1];
rz(-0.22290454) q[1];
rz(1.5570693) q[3];
sx q[3];
rz(-1.1662332) q[3];
sx q[3];
rz(1.3409529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5726418) q[2];
sx q[2];
rz(-0.53223842) q[2];
sx q[2];
rz(0.66298318) q[2];
rz(-0.7524544) q[3];
sx q[3];
rz(-1.3198676) q[3];
sx q[3];
rz(-2.9451356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74333423) q[0];
sx q[0];
rz(-0.4011811) q[0];
sx q[0];
rz(2.4617526) q[0];
rz(2.4196449) q[1];
sx q[1];
rz(-0.55689055) q[1];
sx q[1];
rz(-0.78071761) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3975383) q[0];
sx q[0];
rz(-0.40073943) q[0];
sx q[0];
rz(1.5512054) q[0];
rz(-pi) q[1];
rz(-1.5268097) q[2];
sx q[2];
rz(-0.086217833) q[2];
sx q[2];
rz(1.1607064) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.36759427) q[1];
sx q[1];
rz(-1.9450608) q[1];
sx q[1];
rz(2.4623929) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2188245) q[3];
sx q[3];
rz(-0.84753643) q[3];
sx q[3];
rz(1.7909602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9146933) q[2];
sx q[2];
rz(-1.9136027) q[2];
sx q[2];
rz(2.4719888) q[2];
rz(-0.90034825) q[3];
sx q[3];
rz(-1.0122296) q[3];
sx q[3];
rz(2.3646234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25903073) q[0];
sx q[0];
rz(-0.43864033) q[0];
sx q[0];
rz(0.90748179) q[0];
rz(-2.5471845) q[1];
sx q[1];
rz(-0.92955697) q[1];
sx q[1];
rz(-1.1297191) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4753455) q[0];
sx q[0];
rz(-2.1985558) q[0];
sx q[0];
rz(-2.9055789) q[0];
x q[1];
rz(0.14950651) q[2];
sx q[2];
rz(-2.0975916) q[2];
sx q[2];
rz(3.0164928) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.10953021) q[1];
sx q[1];
rz(-0.9865762) q[1];
sx q[1];
rz(-1.8971127) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3334585) q[3];
sx q[3];
rz(-0.39565797) q[3];
sx q[3];
rz(-1.3612703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0323459) q[2];
sx q[2];
rz(-1.3999908) q[2];
sx q[2];
rz(-2.4692811) q[2];
rz(-0.0191056) q[3];
sx q[3];
rz(-0.34135434) q[3];
sx q[3];
rz(-0.13489558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30447176) q[0];
sx q[0];
rz(-0.79455513) q[0];
sx q[0];
rz(-2.9735612) q[0];
rz(-1.9145603) q[1];
sx q[1];
rz(-1.9214168) q[1];
sx q[1];
rz(-2.8438445) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15140238) q[0];
sx q[0];
rz(-2.08648) q[0];
sx q[0];
rz(-0.97008743) q[0];
rz(-pi) q[1];
rz(-2.1476168) q[2];
sx q[2];
rz(-1.0712396) q[2];
sx q[2];
rz(0.039488878) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7682957) q[1];
sx q[1];
rz(-1.3650286) q[1];
sx q[1];
rz(-1.5453669) q[1];
x q[2];
rz(3.1409114) q[3];
sx q[3];
rz(-1.7295803) q[3];
sx q[3];
rz(0.53631594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1905404) q[2];
sx q[2];
rz(-2.0420044) q[2];
sx q[2];
rz(-0.61868787) q[2];
rz(-2.3330073) q[3];
sx q[3];
rz(-0.4777258) q[3];
sx q[3];
rz(1.3396858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.488778) q[0];
sx q[0];
rz(-1.6365016) q[0];
sx q[0];
rz(-2.6763647) q[0];
rz(-2.6564927) q[1];
sx q[1];
rz(-0.41050375) q[1];
sx q[1];
rz(3.0533275) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7951435) q[0];
sx q[0];
rz(-0.69899054) q[0];
sx q[0];
rz(2.3153852) q[0];
rz(1.5563929) q[2];
sx q[2];
rz(-0.72897899) q[2];
sx q[2];
rz(2.5165503) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8264721) q[1];
sx q[1];
rz(-1.7355647) q[1];
sx q[1];
rz(2.3830448) q[1];
rz(-0.28782423) q[3];
sx q[3];
rz(-1.6555641) q[3];
sx q[3];
rz(0.945795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7469067) q[2];
sx q[2];
rz(-2.2862819) q[2];
sx q[2];
rz(0.49760231) q[2];
rz(-2.6615214) q[3];
sx q[3];
rz(-0.92104715) q[3];
sx q[3];
rz(0.068643071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29480544) q[0];
sx q[0];
rz(-2.6755896) q[0];
sx q[0];
rz(1.71126) q[0];
rz(1.3149186) q[1];
sx q[1];
rz(-0.39350915) q[1];
sx q[1];
rz(1.5865145) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67051149) q[0];
sx q[0];
rz(-1.9626612) q[0];
sx q[0];
rz(3.0125822) q[0];
rz(-pi) q[1];
rz(-2.2261593) q[2];
sx q[2];
rz(-0.829773) q[2];
sx q[2];
rz(-2.1754153) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5615214) q[1];
sx q[1];
rz(-2.5513756) q[1];
sx q[1];
rz(-1.7529704) q[1];
rz(-pi) q[2];
rz(0.076818941) q[3];
sx q[3];
rz(-0.94074501) q[3];
sx q[3];
rz(-2.1253642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2985208) q[2];
sx q[2];
rz(-2.4303747) q[2];
sx q[2];
rz(-2.2141875) q[2];
rz(-1.1585506) q[3];
sx q[3];
rz(-0.68803334) q[3];
sx q[3];
rz(-3.0454175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7189099) q[0];
sx q[0];
rz(-2.2490608) q[0];
sx q[0];
rz(-2.6655777) q[0];
rz(-2.1503275) q[1];
sx q[1];
rz(-2.3920993) q[1];
sx q[1];
rz(3.057726) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80591969) q[0];
sx q[0];
rz(-1.085338) q[0];
sx q[0];
rz(0.40747633) q[0];
rz(1.9344566) q[2];
sx q[2];
rz(-1.1354828) q[2];
sx q[2];
rz(-0.93863034) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8955651) q[1];
sx q[1];
rz(-1.9894587) q[1];
sx q[1];
rz(-3.0484867) q[1];
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
rz(2.2241761) q[2];
sx q[2];
rz(-0.26599628) q[2];
sx q[2];
rz(0.64845294) q[2];
rz(-2.2367541) q[3];
sx q[3];
rz(-1.4584533) q[3];
sx q[3];
rz(-2.8308433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34173486) q[0];
sx q[0];
rz(-0.71001995) q[0];
sx q[0];
rz(-0.56890666) q[0];
rz(2.3078602) q[1];
sx q[1];
rz(-1.8582452) q[1];
sx q[1];
rz(-1.0047097) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89812714) q[0];
sx q[0];
rz(-1.3782189) q[0];
sx q[0];
rz(0.63611998) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23022454) q[2];
sx q[2];
rz(-2.0966999) q[2];
sx q[2];
rz(0.81947015) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8160287) q[1];
sx q[1];
rz(-2.4392011) q[1];
sx q[1];
rz(1.4588474) q[1];
rz(-pi) q[2];
rz(-2.5540657) q[3];
sx q[3];
rz(-1.5024795) q[3];
sx q[3];
rz(2.5053566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7052762) q[2];
sx q[2];
rz(-2.269561) q[2];
sx q[2];
rz(-0.3479859) q[2];
rz(2.3157388) q[3];
sx q[3];
rz(-2.9233942) q[3];
sx q[3];
rz(2.5364449) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0843049) q[0];
sx q[0];
rz(-1.0658406) q[0];
sx q[0];
rz(2.9344946) q[0];
rz(-0.53889489) q[1];
sx q[1];
rz(-0.62108827) q[1];
sx q[1];
rz(-2.5394687) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83610632) q[0];
sx q[0];
rz(-1.0998816) q[0];
sx q[0];
rz(0.54022809) q[0];
rz(1.978157) q[2];
sx q[2];
rz(-1.0216139) q[2];
sx q[2];
rz(-2.801535) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.45710292) q[1];
sx q[1];
rz(-2.7278363) q[1];
sx q[1];
rz(1.7587508) q[1];
rz(1.8982417) q[3];
sx q[3];
rz(-2.2909938) q[3];
sx q[3];
rz(2.7011288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.78486097) q[2];
sx q[2];
rz(-2.1897903) q[2];
sx q[2];
rz(-1.7968563) q[2];
rz(2.6209659) q[3];
sx q[3];
rz(-2.8713363) q[3];
sx q[3];
rz(0.80210137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.44705924) q[0];
sx q[0];
rz(-2.4166528) q[0];
sx q[0];
rz(1.7550533) q[0];
rz(-2.3583892) q[1];
sx q[1];
rz(-1.422311) q[1];
sx q[1];
rz(1.5625988) q[1];
rz(-1.8405909) q[2];
sx q[2];
rz(-1.2392344) q[2];
sx q[2];
rz(0.45891529) q[2];
rz(-0.047719638) q[3];
sx q[3];
rz(-0.24038355) q[3];
sx q[3];
rz(3.0265831) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
