OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1908258) q[0];
sx q[0];
rz(-1.289225) q[0];
sx q[0];
rz(3.0486795) q[0];
rz(3.0463123) q[1];
sx q[1];
rz(-2.4089101) q[1];
sx q[1];
rz(-1.9021775) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5372972) q[0];
sx q[0];
rz(-1.3193466) q[0];
sx q[0];
rz(2.9049113) q[0];
rz(-pi) q[1];
rz(-0.26387604) q[2];
sx q[2];
rz(-1.7603163) q[2];
sx q[2];
rz(-2.5352728) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.79151407) q[1];
sx q[1];
rz(-2.8100612) q[1];
sx q[1];
rz(-0.091011957) q[1];
rz(0.69783437) q[3];
sx q[3];
rz(-2.4437332) q[3];
sx q[3];
rz(2.7891911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6713082) q[2];
sx q[2];
rz(-1.6613864) q[2];
sx q[2];
rz(1.3489464) q[2];
rz(-0.74364439) q[3];
sx q[3];
rz(-2.8642004) q[3];
sx q[3];
rz(-0.17717895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.141356) q[0];
sx q[0];
rz(-1.5445222) q[0];
sx q[0];
rz(-0.29362383) q[0];
rz(-2.1108421) q[1];
sx q[1];
rz(-1.7799957) q[1];
sx q[1];
rz(-0.34293276) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6688419) q[0];
sx q[0];
rz(-1.051051) q[0];
sx q[0];
rz(-0.23730554) q[0];
rz(-pi) q[1];
rz(0.12983506) q[2];
sx q[2];
rz(-1.9399683) q[2];
sx q[2];
rz(-0.8952199) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1520815) q[1];
sx q[1];
rz(-0.73484269) q[1];
sx q[1];
rz(0.32232018) q[1];
x q[2];
rz(2.8906451) q[3];
sx q[3];
rz(-2.1852369) q[3];
sx q[3];
rz(-2.3411103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.12884101) q[2];
sx q[2];
rz(-1.7495456) q[2];
sx q[2];
rz(-0.7592321) q[2];
rz(3.1208842) q[3];
sx q[3];
rz(-1.2660675) q[3];
sx q[3];
rz(0.43829632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.6193806) q[0];
sx q[0];
rz(-2.5396357) q[0];
sx q[0];
rz(0.0044599175) q[0];
rz(0.64747539) q[1];
sx q[1];
rz(-2.5471893) q[1];
sx q[1];
rz(-0.78537816) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4606708) q[0];
sx q[0];
rz(-2.9630532) q[0];
sx q[0];
rz(2.6723249) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0702281) q[2];
sx q[2];
rz(-2.7874261) q[2];
sx q[2];
rz(-0.11812299) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8632311) q[1];
sx q[1];
rz(-1.0309217) q[1];
sx q[1];
rz(-3.0616772) q[1];
rz(-pi) q[2];
x q[2];
rz(2.249751) q[3];
sx q[3];
rz(-1.0247318) q[3];
sx q[3];
rz(-0.31072703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.37041) q[2];
sx q[2];
rz(-2.2213171) q[2];
sx q[2];
rz(1.3578337) q[2];
rz(-2.8010098) q[3];
sx q[3];
rz(-2.1850696) q[3];
sx q[3];
rz(-0.79184872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99712813) q[0];
sx q[0];
rz(-2.0576532) q[0];
sx q[0];
rz(-0.88515627) q[0];
rz(-0.74686933) q[1];
sx q[1];
rz(-0.62364686) q[1];
sx q[1];
rz(2.0451827) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88056662) q[0];
sx q[0];
rz(-1.2231022) q[0];
sx q[0];
rz(2.8096922) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0463767) q[2];
sx q[2];
rz(-0.066844373) q[2];
sx q[2];
rz(-2.9796556) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4064286) q[1];
sx q[1];
rz(-0.53710912) q[1];
sx q[1];
rz(0.5237315) q[1];
rz(2.0799594) q[3];
sx q[3];
rz(-2.235755) q[3];
sx q[3];
rz(3.0388415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5235644) q[2];
sx q[2];
rz(-2.9310493) q[2];
sx q[2];
rz(-3.1169685) q[2];
rz(2.5060182) q[3];
sx q[3];
rz(-0.98819757) q[3];
sx q[3];
rz(-0.88328254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4516975) q[0];
sx q[0];
rz(-2.8391333) q[0];
sx q[0];
rz(-2.6746993) q[0];
rz(-3.0310071) q[1];
sx q[1];
rz(-2.6778335) q[1];
sx q[1];
rz(-2.0707524) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3239646) q[0];
sx q[0];
rz(-1.1187129) q[0];
sx q[0];
rz(0.26497193) q[0];
rz(-pi) q[1];
rz(-2.8230225) q[2];
sx q[2];
rz(-1.7442707) q[2];
sx q[2];
rz(-1.3248548) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5731949) q[1];
sx q[1];
rz(-1.0612951) q[1];
sx q[1];
rz(-1.1361109) q[1];
rz(-0.91584622) q[3];
sx q[3];
rz(-0.98036843) q[3];
sx q[3];
rz(-0.46186033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0232627) q[2];
sx q[2];
rz(-1.3619962) q[2];
sx q[2];
rz(2.2892717) q[2];
rz(0.037633745) q[3];
sx q[3];
rz(-1.9084385) q[3];
sx q[3];
rz(-2.5467303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1325876) q[0];
sx q[0];
rz(-1.9629033) q[0];
sx q[0];
rz(1.8527385) q[0];
rz(1.6291078) q[1];
sx q[1];
rz(-1.1237203) q[1];
sx q[1];
rz(1.1423133) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70002551) q[0];
sx q[0];
rz(-1.6684932) q[0];
sx q[0];
rz(-2.7872681) q[0];
rz(-pi) q[1];
x q[1];
rz(1.283848) q[2];
sx q[2];
rz(-2.9252238) q[2];
sx q[2];
rz(-1.8277825) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9276322) q[1];
sx q[1];
rz(-2.7134656) q[1];
sx q[1];
rz(0.63251782) q[1];
x q[2];
rz(-0.31123881) q[3];
sx q[3];
rz(-1.7390307) q[3];
sx q[3];
rz(2.5995035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3219354) q[2];
sx q[2];
rz(-1.8737996) q[2];
sx q[2];
rz(-3.0211871) q[2];
rz(-1.985792) q[3];
sx q[3];
rz(-1.6121696) q[3];
sx q[3];
rz(-0.22404484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33893809) q[0];
sx q[0];
rz(-1.3405565) q[0];
sx q[0];
rz(1.6876203) q[0];
rz(1.5441719) q[1];
sx q[1];
rz(-1.6463966) q[1];
sx q[1];
rz(-2.6905751) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9748443) q[0];
sx q[0];
rz(-1.6303501) q[0];
sx q[0];
rz(2.8400499) q[0];
x q[1];
rz(-1.9372378) q[2];
sx q[2];
rz(-1.7882573) q[2];
sx q[2];
rz(1.5254453) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9684194) q[1];
sx q[1];
rz(-2.0035158) q[1];
sx q[1];
rz(1.7041901) q[1];
rz(-pi) q[2];
rz(-2.3312206) q[3];
sx q[3];
rz(-1.5640508) q[3];
sx q[3];
rz(-0.41616671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1958127) q[2];
sx q[2];
rz(-1.7273629) q[2];
sx q[2];
rz(1.3607402) q[2];
rz(-2.1360548) q[3];
sx q[3];
rz(-1.9660299) q[3];
sx q[3];
rz(1.7520693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1239531) q[0];
sx q[0];
rz(-3.0160976) q[0];
sx q[0];
rz(2.4355167) q[0];
rz(1.5438682) q[1];
sx q[1];
rz(-1.7192625) q[1];
sx q[1];
rz(-0.85618883) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5917011) q[0];
sx q[0];
rz(-0.020792637) q[0];
sx q[0];
rz(-1.2077622) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43120877) q[2];
sx q[2];
rz(-0.87006205) q[2];
sx q[2];
rz(1.7675811) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.013102839) q[1];
sx q[1];
rz(-2.1399413) q[1];
sx q[1];
rz(-1.3291969) q[1];
rz(-1.8185577) q[3];
sx q[3];
rz(-2.6399586) q[3];
sx q[3];
rz(0.35546965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2989444) q[2];
sx q[2];
rz(-1.4010669) q[2];
sx q[2];
rz(-2.974158) q[2];
rz(1.5542479) q[3];
sx q[3];
rz(-2.3685679) q[3];
sx q[3];
rz(-2.1412444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-1.7635968) q[0];
sx q[0];
rz(-1.4507699) q[0];
sx q[0];
rz(0.3717306) q[0];
rz(-1.7077712) q[1];
sx q[1];
rz(-1.5377518) q[1];
sx q[1];
rz(-0.30002123) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8512631) q[0];
sx q[0];
rz(-1.6870878) q[0];
sx q[0];
rz(2.854748) q[0];
rz(-pi) q[1];
rz(-0.062252684) q[2];
sx q[2];
rz(-0.57800284) q[2];
sx q[2];
rz(-2.5684772) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.91050657) q[1];
sx q[1];
rz(-2.4795737) q[1];
sx q[1];
rz(2.9390017) q[1];
x q[2];
rz(-2.2737502) q[3];
sx q[3];
rz(-1.6658837) q[3];
sx q[3];
rz(0.33012182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6649449) q[2];
sx q[2];
rz(-2.6740394) q[2];
sx q[2];
rz(0.42759582) q[2];
rz(-1.5852196) q[3];
sx q[3];
rz(-2.0245602) q[3];
sx q[3];
rz(-2.3868886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72117358) q[0];
sx q[0];
rz(-0.54062802) q[0];
sx q[0];
rz(2.6721201) q[0];
rz(2.2162614) q[1];
sx q[1];
rz(-2.053849) q[1];
sx q[1];
rz(-3.1184149) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3566475) q[0];
sx q[0];
rz(-2.0933505) q[0];
sx q[0];
rz(1.7607259) q[0];
rz(-pi) q[1];
rz(2.542664) q[2];
sx q[2];
rz(-0.99236503) q[2];
sx q[2];
rz(-2.3430062) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7851023) q[1];
sx q[1];
rz(-1.7447326) q[1];
sx q[1];
rz(-0.16696232) q[1];
rz(-2.9523207) q[3];
sx q[3];
rz(-0.84976518) q[3];
sx q[3];
rz(1.4355575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.8699441) q[2];
sx q[2];
rz(-1.2056489) q[2];
sx q[2];
rz(0.49087697) q[2];
rz(2.0513746) q[3];
sx q[3];
rz(-2.1722983) q[3];
sx q[3];
rz(0.31392613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8563817) q[0];
sx q[0];
rz(-0.87775341) q[0];
sx q[0];
rz(-2.0650771) q[0];
rz(0.27676997) q[1];
sx q[1];
rz(-2.3634187) q[1];
sx q[1];
rz(1.707911) q[1];
rz(-1.8770915) q[2];
sx q[2];
rz(-2.2363792) q[2];
sx q[2];
rz(-2.1909942) q[2];
rz(-3.1076486) q[3];
sx q[3];
rz(-2.6068519) q[3];
sx q[3];
rz(3.1153771) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
