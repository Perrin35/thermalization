OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.45286173) q[0];
sx q[0];
rz(2.9714669) q[0];
sx q[0];
rz(10.210769) q[0];
rz(0.6056447) q[1];
sx q[1];
rz(-2.4906467) q[1];
sx q[1];
rz(-0.63408607) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68361002) q[0];
sx q[0];
rz(-1.5994706) q[0];
sx q[0];
rz(0.81575127) q[0];
x q[1];
rz(1.6526821) q[2];
sx q[2];
rz(-2.4856644) q[2];
sx q[2];
rz(-1.3002849) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4418728) q[1];
sx q[1];
rz(-0.77077121) q[1];
sx q[1];
rz(2.2295879) q[1];
rz(-pi) q[2];
rz(0.4001873) q[3];
sx q[3];
rz(-1.549198) q[3];
sx q[3];
rz(0.17822972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2259851) q[2];
sx q[2];
rz(-2.6245124) q[2];
sx q[2];
rz(-1.8784286) q[2];
rz(-1.2849215) q[3];
sx q[3];
rz(-1.4611171) q[3];
sx q[3];
rz(-0.29299709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15853515) q[0];
sx q[0];
rz(-1.2936658) q[0];
sx q[0];
rz(-2.7040226) q[0];
rz(0.63105398) q[1];
sx q[1];
rz(-0.32866207) q[1];
sx q[1];
rz(0.22110573) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3723345) q[0];
sx q[0];
rz(-0.49494574) q[0];
sx q[0];
rz(-1.9702205) q[0];
x q[1];
rz(-0.11110335) q[2];
sx q[2];
rz(-1.2239893) q[2];
sx q[2];
rz(2.8135516) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.088614956) q[1];
sx q[1];
rz(-0.47514519) q[1];
sx q[1];
rz(-0.26762025) q[1];
rz(-pi) q[2];
rz(-2.9651871) q[3];
sx q[3];
rz(-0.50900148) q[3];
sx q[3];
rz(-2.727946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.21800403) q[2];
sx q[2];
rz(-1.4322586) q[2];
sx q[2];
rz(2.5533) q[2];
rz(-0.44899392) q[3];
sx q[3];
rz(-0.42789999) q[3];
sx q[3];
rz(2.0935521) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5730729) q[0];
sx q[0];
rz(-3.0439215) q[0];
sx q[0];
rz(2.1411238) q[0];
rz(0.72552848) q[1];
sx q[1];
rz(-1.0745445) q[1];
sx q[1];
rz(-0.75769889) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.476186) q[0];
sx q[0];
rz(-1.889567) q[0];
sx q[0];
rz(-0.55253367) q[0];
x q[1];
rz(-1.7311086) q[2];
sx q[2];
rz(-1.7756887) q[2];
sx q[2];
rz(0.93970539) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.11848395) q[1];
sx q[1];
rz(-0.99220905) q[1];
sx q[1];
rz(0.56225496) q[1];
x q[2];
rz(-2.9919639) q[3];
sx q[3];
rz(-0.64390874) q[3];
sx q[3];
rz(0.097188918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9329325) q[2];
sx q[2];
rz(-2.9186086) q[2];
sx q[2];
rz(-0.83646742) q[2];
rz(1.4423192) q[3];
sx q[3];
rz(-1.6675555) q[3];
sx q[3];
rz(-2.9377655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1214685) q[0];
sx q[0];
rz(-3.0614873) q[0];
sx q[0];
rz(-2.6304723) q[0];
rz(-3.064149) q[1];
sx q[1];
rz(-0.42235342) q[1];
sx q[1];
rz(-1.8130594) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9851345) q[0];
sx q[0];
rz(-1.7256323) q[0];
sx q[0];
rz(1.3550718) q[0];
rz(0.76339108) q[2];
sx q[2];
rz(-1.1995458) q[2];
sx q[2];
rz(3.0766746) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7718539) q[1];
sx q[1];
rz(-0.26680294) q[1];
sx q[1];
rz(1.5727444) q[1];
x q[2];
rz(-1.3489181) q[3];
sx q[3];
rz(-1.4639877) q[3];
sx q[3];
rz(0.91339236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0478583) q[2];
sx q[2];
rz(-1.2183943) q[2];
sx q[2];
rz(1.1070586) q[2];
rz(-0.47248653) q[3];
sx q[3];
rz(-1.6146336) q[3];
sx q[3];
rz(-2.2842177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1012786) q[0];
sx q[0];
rz(-0.89656985) q[0];
sx q[0];
rz(2.8033946) q[0];
rz(1.8473373) q[1];
sx q[1];
rz(-1.5599524) q[1];
sx q[1];
rz(2.6370874) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38406661) q[0];
sx q[0];
rz(-2.5528918) q[0];
sx q[0];
rz(-1.9996044) q[0];
rz(-pi) q[1];
rz(-1.9485103) q[2];
sx q[2];
rz(-2.0078805) q[2];
sx q[2];
rz(-1.9405685) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.63197631) q[1];
sx q[1];
rz(-2.7512433) q[1];
sx q[1];
rz(0.901464) q[1];
rz(-pi) q[2];
rz(-0.65977804) q[3];
sx q[3];
rz(-2.1187083) q[3];
sx q[3];
rz(-1.1783311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.111104) q[2];
sx q[2];
rz(-0.8384853) q[2];
sx q[2];
rz(1.6476691) q[2];
rz(1.4533639) q[3];
sx q[3];
rz(-2.2036392) q[3];
sx q[3];
rz(1.7593613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9313653) q[0];
sx q[0];
rz(-2.2670822) q[0];
sx q[0];
rz(2.0507623) q[0];
rz(0.53972721) q[1];
sx q[1];
rz(-1.9878186) q[1];
sx q[1];
rz(2.9528023) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8180346) q[0];
sx q[0];
rz(-1.8206882) q[0];
sx q[0];
rz(-0.37139335) q[0];
x q[1];
rz(1.9300869) q[2];
sx q[2];
rz(-0.28887666) q[2];
sx q[2];
rz(-1.2668244) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.69953883) q[1];
sx q[1];
rz(-0.25146723) q[1];
sx q[1];
rz(0.26775189) q[1];
x q[2];
rz(-0.24174989) q[3];
sx q[3];
rz(-1.3931735) q[3];
sx q[3];
rz(-2.1135981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9699036) q[2];
sx q[2];
rz(-2.1128555) q[2];
sx q[2];
rz(1.3151273) q[2];
rz(1.5054437) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(1.858254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0514907) q[0];
sx q[0];
rz(-0.83054709) q[0];
sx q[0];
rz(0.32399696) q[0];
rz(-1.3011159) q[1];
sx q[1];
rz(-2.3062861) q[1];
sx q[1];
rz(1.3791929) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8934879) q[0];
sx q[0];
rz(-1.6525434) q[0];
sx q[0];
rz(-1.8155314) q[0];
x q[1];
rz(-0.46160134) q[2];
sx q[2];
rz(-1.201655) q[2];
sx q[2];
rz(-2.387407) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0566237) q[1];
sx q[1];
rz(-1.8693722) q[1];
sx q[1];
rz(1.741239) q[1];
rz(-pi) q[2];
rz(-2.561065) q[3];
sx q[3];
rz(-1.0167828) q[3];
sx q[3];
rz(1.1298657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1088915) q[2];
sx q[2];
rz(-1.6872493) q[2];
sx q[2];
rz(0.94318715) q[2];
rz(0.33106783) q[3];
sx q[3];
rz(-1.7739242) q[3];
sx q[3];
rz(-2.846472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11809764) q[0];
sx q[0];
rz(-0.58693111) q[0];
sx q[0];
rz(-0.76422894) q[0];
rz(3.0006192) q[1];
sx q[1];
rz(-0.39607513) q[1];
sx q[1];
rz(-1.2932628) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30771502) q[0];
sx q[0];
rz(-2.6617962) q[0];
sx q[0];
rz(-1.1307554) q[0];
rz(-pi) q[1];
rz(0.60111945) q[2];
sx q[2];
rz(-1.4919315) q[2];
sx q[2];
rz(-2.2476946) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5544719) q[1];
sx q[1];
rz(-1.9738102) q[1];
sx q[1];
rz(2.9906669) q[1];
rz(1.3702277) q[3];
sx q[3];
rz(-1.1784394) q[3];
sx q[3];
rz(1.4810824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76688898) q[2];
sx q[2];
rz(-2.1024487) q[2];
sx q[2];
rz(2.5602706) q[2];
rz(-2.2733722) q[3];
sx q[3];
rz(-0.41728443) q[3];
sx q[3];
rz(-1.8301615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54810369) q[0];
sx q[0];
rz(-0.5797525) q[0];
sx q[0];
rz(1.2506437) q[0];
rz(-0.99682322) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(1.8539799) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33289136) q[0];
sx q[0];
rz(-1.5960072) q[0];
sx q[0];
rz(2.0602134) q[0];
rz(-1.397923) q[2];
sx q[2];
rz(-1.0908974) q[2];
sx q[2];
rz(-2.929504) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1907562) q[1];
sx q[1];
rz(-1.7777998) q[1];
sx q[1];
rz(-0.91916577) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0887293) q[3];
sx q[3];
rz(-0.82517805) q[3];
sx q[3];
rz(2.4406274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.320497) q[2];
sx q[2];
rz(-0.38791502) q[2];
sx q[2];
rz(-1.256475) q[2];
rz(-0.16658941) q[3];
sx q[3];
rz(-1.5649256) q[3];
sx q[3];
rz(-1.0725718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88084108) q[0];
sx q[0];
rz(-2.3374225) q[0];
sx q[0];
rz(2.8163731) q[0];
rz(1.1351769) q[1];
sx q[1];
rz(-2.4644641) q[1];
sx q[1];
rz(0.36718711) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67883915) q[0];
sx q[0];
rz(-0.05719962) q[0];
sx q[0];
rz(0.57533933) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2904097) q[2];
sx q[2];
rz(-2.1791611) q[2];
sx q[2];
rz(2.8031363) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2314184) q[1];
sx q[1];
rz(-1.7703729) q[1];
sx q[1];
rz(2.477596) q[1];
rz(-0.48093421) q[3];
sx q[3];
rz(-2.1554865) q[3];
sx q[3];
rz(1.9422873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8964768) q[2];
sx q[2];
rz(-0.8090691) q[2];
sx q[2];
rz(2.0754576) q[2];
rz(0.079244763) q[3];
sx q[3];
rz(-0.60278046) q[3];
sx q[3];
rz(1.4762896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8463678) q[0];
sx q[0];
rz(-1.5548779) q[0];
sx q[0];
rz(1.5665733) q[0];
rz(0.29905839) q[1];
sx q[1];
rz(-2.538264) q[1];
sx q[1];
rz(-2.6187142) q[1];
rz(0.085061442) q[2];
sx q[2];
rz(-0.75123514) q[2];
sx q[2];
rz(-3.0820465) q[2];
rz(-0.15300898) q[3];
sx q[3];
rz(-2.2867793) q[3];
sx q[3];
rz(2.2147562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];