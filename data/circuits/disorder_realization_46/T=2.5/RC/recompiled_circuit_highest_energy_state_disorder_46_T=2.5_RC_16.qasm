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
rz(-2.6645339) q[0];
sx q[0];
rz(-2.8395489) q[0];
sx q[0];
rz(-1.8855236) q[0];
rz(0.37580252) q[1];
sx q[1];
rz(1.8355651) q[1];
sx q[1];
rz(10.553283) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1953729) q[0];
sx q[0];
rz(-1.8502619) q[0];
sx q[0];
rz(-2.165876) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8774421) q[2];
sx q[2];
rz(-2.1696343) q[2];
sx q[2];
rz(-1.6055589) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.30762954) q[1];
sx q[1];
rz(-1.071573) q[1];
sx q[1];
rz(-2.7021033) q[1];
rz(2.5735568) q[3];
sx q[3];
rz(-2.4246693) q[3];
sx q[3];
rz(2.4942945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6065373) q[2];
sx q[2];
rz(-2.8474319) q[2];
sx q[2];
rz(1.199022) q[2];
rz(2.8626275) q[3];
sx q[3];
rz(-1.4594892) q[3];
sx q[3];
rz(0.038486686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26175427) q[0];
sx q[0];
rz(-0.45670515) q[0];
sx q[0];
rz(0.36619151) q[0];
rz(1.2884033) q[1];
sx q[1];
rz(-2.2346456) q[1];
sx q[1];
rz(0.2624951) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27124559) q[0];
sx q[0];
rz(-2.7923005) q[0];
sx q[0];
rz(0.55089921) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0832564) q[2];
sx q[2];
rz(-1.3223786) q[2];
sx q[2];
rz(0.89983856) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9001635) q[1];
sx q[1];
rz(-1.2821322) q[1];
sx q[1];
rz(-2.7735387) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6672584) q[3];
sx q[3];
rz(-0.25611195) q[3];
sx q[3];
rz(-0.90463582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.75258201) q[2];
sx q[2];
rz(-1.8730619) q[2];
sx q[2];
rz(-2.9627723) q[2];
rz(-3.1255152) q[3];
sx q[3];
rz(-0.5355081) q[3];
sx q[3];
rz(0.9308365) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43612424) q[0];
sx q[0];
rz(-3.0146764) q[0];
sx q[0];
rz(-0.57417589) q[0];
rz(-0.15159675) q[1];
sx q[1];
rz(-2.6972289) q[1];
sx q[1];
rz(1.9134329) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3686217) q[0];
sx q[0];
rz(-0.54062343) q[0];
sx q[0];
rz(1.878488) q[0];
rz(-1.6723694) q[2];
sx q[2];
rz(-2.5463062) q[2];
sx q[2];
rz(3.1208619) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2252402) q[1];
sx q[1];
rz(-0.95969363) q[1];
sx q[1];
rz(3.0322187) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0869998) q[3];
sx q[3];
rz(-1.3143149) q[3];
sx q[3];
rz(-0.90638049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.14153081) q[2];
sx q[2];
rz(-2.5059293) q[2];
sx q[2];
rz(-0.89111152) q[2];
rz(-2.7530503) q[3];
sx q[3];
rz(-1.5107692) q[3];
sx q[3];
rz(-2.8019606) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5163088) q[0];
sx q[0];
rz(-3.0966274) q[0];
sx q[0];
rz(0.48211023) q[0];
rz(2.4730143) q[1];
sx q[1];
rz(-1.6079638) q[1];
sx q[1];
rz(-0.63749981) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9264049) q[0];
sx q[0];
rz(-0.63894659) q[0];
sx q[0];
rz(1.9685345) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0257865) q[2];
sx q[2];
rz(-2.2409332) q[2];
sx q[2];
rz(1.1095123) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1177534) q[1];
sx q[1];
rz(-1.6477403) q[1];
sx q[1];
rz(1.5368098) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25084514) q[3];
sx q[3];
rz(-1.224784) q[3];
sx q[3];
rz(0.097867997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.94347) q[2];
sx q[2];
rz(-2.1169457) q[2];
sx q[2];
rz(-0.060297273) q[2];
rz(0.30334011) q[3];
sx q[3];
rz(-0.077849418) q[3];
sx q[3];
rz(0.2651324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9918793) q[0];
sx q[0];
rz(-2.7880221) q[0];
sx q[0];
rz(0.41325945) q[0];
rz(2.7464271) q[1];
sx q[1];
rz(-1.7839849) q[1];
sx q[1];
rz(-1.0023592) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7089139) q[0];
sx q[0];
rz(-1.029338) q[0];
sx q[0];
rz(1.0546636) q[0];
x q[1];
rz(-0.036040739) q[2];
sx q[2];
rz(-1.7805575) q[2];
sx q[2];
rz(-2.710611) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.614568) q[1];
sx q[1];
rz(-0.90674149) q[1];
sx q[1];
rz(-0.27399691) q[1];
rz(-pi) q[2];
rz(-1.0071383) q[3];
sx q[3];
rz(-0.37630577) q[3];
sx q[3];
rz(0.2175771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5137382) q[2];
sx q[2];
rz(-2.0444874) q[2];
sx q[2];
rz(1.537079) q[2];
rz(-1.1995992) q[3];
sx q[3];
rz(-1.1211841) q[3];
sx q[3];
rz(-0.77170038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5143249) q[0];
sx q[0];
rz(-1.6384614) q[0];
sx q[0];
rz(-2.7728873) q[0];
rz(-0.74327028) q[1];
sx q[1];
rz(-0.5629881) q[1];
sx q[1];
rz(-0.35400131) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0403772) q[0];
sx q[0];
rz(-0.64283744) q[0];
sx q[0];
rz(1.6580216) q[0];
x q[1];
rz(-1.8839624) q[2];
sx q[2];
rz(-1.0651565) q[2];
sx q[2];
rz(-2.2693116) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.88601513) q[1];
sx q[1];
rz(-1.3929345) q[1];
sx q[1];
rz(-2.1484003) q[1];
rz(0.50320585) q[3];
sx q[3];
rz(-0.70069289) q[3];
sx q[3];
rz(-1.038674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1165498) q[2];
sx q[2];
rz(-1.9321238) q[2];
sx q[2];
rz(-3.1361191) q[2];
rz(-2.6239851) q[3];
sx q[3];
rz(-0.36492437) q[3];
sx q[3];
rz(0.027211729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40590498) q[0];
sx q[0];
rz(-2.5820177) q[0];
sx q[0];
rz(0.21301183) q[0];
rz(1.4951911) q[1];
sx q[1];
rz(-2.5081684) q[1];
sx q[1];
rz(2.2921553) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0606468) q[0];
sx q[0];
rz(-2.1889926) q[0];
sx q[0];
rz(-1.5175876) q[0];
rz(0.030730129) q[2];
sx q[2];
rz(-1.5275914) q[2];
sx q[2];
rz(1.876056) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.53687261) q[1];
sx q[1];
rz(-1.1013796) q[1];
sx q[1];
rz(-2.0783271) q[1];
rz(-pi) q[2];
rz(-1.0424741) q[3];
sx q[3];
rz(-1.436621) q[3];
sx q[3];
rz(-3.1229916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9269632) q[2];
sx q[2];
rz(-1.9081076) q[2];
sx q[2];
rz(0.54369533) q[2];
rz(-2.8804273) q[3];
sx q[3];
rz(-1.5615014) q[3];
sx q[3];
rz(-0.15682596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1866622) q[0];
sx q[0];
rz(-1.0533227) q[0];
sx q[0];
rz(0.18761158) q[0];
rz(-1.1232417) q[1];
sx q[1];
rz(-0.76485991) q[1];
sx q[1];
rz(0.16256464) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.150972) q[0];
sx q[0];
rz(-0.91998581) q[0];
sx q[0];
rz(1.1605754) q[0];
rz(-1.8674054) q[2];
sx q[2];
rz(-0.54365082) q[2];
sx q[2];
rz(0.44270502) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.39531884) q[1];
sx q[1];
rz(-0.95764179) q[1];
sx q[1];
rz(-0.055387251) q[1];
rz(-pi) q[2];
rz(-3.0288234) q[3];
sx q[3];
rz(-0.9774607) q[3];
sx q[3];
rz(1.6360374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2742013) q[2];
sx q[2];
rz(-2.826639) q[2];
sx q[2];
rz(-2.9401722) q[2];
rz(-2.6650186) q[3];
sx q[3];
rz(-1.7628935) q[3];
sx q[3];
rz(-0.99982888) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1808712) q[0];
sx q[0];
rz(-0.23715401) q[0];
sx q[0];
rz(0.57688212) q[0];
rz(-0.30413973) q[1];
sx q[1];
rz(-1.1241333) q[1];
sx q[1];
rz(-1.1041799) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3873515) q[0];
sx q[0];
rz(-1.9490483) q[0];
sx q[0];
rz(-3.0463112) q[0];
rz(-pi) q[1];
rz(2.9530802) q[2];
sx q[2];
rz(-2.5645263) q[2];
sx q[2];
rz(-0.12332502) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6526281) q[1];
sx q[1];
rz(-1.1833937) q[1];
sx q[1];
rz(2.1968958) q[1];
rz(-pi) q[2];
rz(-1.1796168) q[3];
sx q[3];
rz(-1.5857539) q[3];
sx q[3];
rz(-2.7103031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1748109) q[2];
sx q[2];
rz(-0.2468144) q[2];
sx q[2];
rz(0.68320572) q[2];
rz(-1.5934058) q[3];
sx q[3];
rz(-0.67845172) q[3];
sx q[3];
rz(-2.9233176) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9504647) q[0];
sx q[0];
rz(-0.8466962) q[0];
sx q[0];
rz(-0.51093131) q[0];
rz(1.9966985) q[1];
sx q[1];
rz(-1.9547918) q[1];
sx q[1];
rz(1.5085545) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.78054) q[0];
sx q[0];
rz(-2.1117438) q[0];
sx q[0];
rz(-1.279328) q[0];
rz(1.0912446) q[2];
sx q[2];
rz(-2.6065718) q[2];
sx q[2];
rz(2.1370691) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3737693) q[1];
sx q[1];
rz(-1.1255742) q[1];
sx q[1];
rz(-1.3740463) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4436117) q[3];
sx q[3];
rz(-2.7578286) q[3];
sx q[3];
rz(-2.021054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7213584) q[2];
sx q[2];
rz(-0.63174641) q[2];
sx q[2];
rz(1.4314123) q[2];
rz(0.015627705) q[3];
sx q[3];
rz(-1.8764007) q[3];
sx q[3];
rz(-2.6354852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(2.2331727) q[0];
sx q[0];
rz(-1.7411727) q[0];
sx q[0];
rz(2.2199051) q[0];
rz(-1.635101) q[1];
sx q[1];
rz(-1.2163305) q[1];
sx q[1];
rz(-0.46179927) q[1];
rz(-0.93585451) q[2];
sx q[2];
rz(-2.6781185) q[2];
sx q[2];
rz(0.92739633) q[2];
rz(1.6484281) q[3];
sx q[3];
rz(-1.6059198) q[3];
sx q[3];
rz(-3.0618659) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
