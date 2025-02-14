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
rz(-1.0128385) q[0];
sx q[0];
rz(-1.9586118) q[0];
sx q[0];
rz(2.8943789) q[0];
rz(1.5565058) q[1];
sx q[1];
rz(-1.0143919) q[1];
sx q[1];
rz(2.1870764) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.010095) q[0];
sx q[0];
rz(-2.7662377) q[0];
sx q[0];
rz(-1.6502871) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45329161) q[2];
sx q[2];
rz(-1.5201836) q[2];
sx q[2];
rz(-2.5091189) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.40084154) q[1];
sx q[1];
rz(-2.2993326) q[1];
sx q[1];
rz(0.42899112) q[1];
rz(-2.8273647) q[3];
sx q[3];
rz(-0.79821051) q[3];
sx q[3];
rz(1.8857005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1859833) q[2];
sx q[2];
rz(-0.60638967) q[2];
sx q[2];
rz(2.8601904) q[2];
rz(1.227281) q[3];
sx q[3];
rz(-1.3325007) q[3];
sx q[3];
rz(-1.9277771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94673741) q[0];
sx q[0];
rz(-0.40732107) q[0];
sx q[0];
rz(-2.5545004) q[0];
rz(-2.6214444) q[1];
sx q[1];
rz(-1.2112434) q[1];
sx q[1];
rz(2.928226) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0555455) q[0];
sx q[0];
rz(-0.16250691) q[0];
sx q[0];
rz(-1.1799728) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93307067) q[2];
sx q[2];
rz(-2.3847849) q[2];
sx q[2];
rz(-0.30100664) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0679464) q[1];
sx q[1];
rz(-1.6151307) q[1];
sx q[1];
rz(1.4774051) q[1];
rz(-2.0265686) q[3];
sx q[3];
rz(-2.5682862) q[3];
sx q[3];
rz(0.70070964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0352036) q[2];
sx q[2];
rz(-0.28855244) q[2];
sx q[2];
rz(2.8698548) q[2];
rz(2.4954097) q[3];
sx q[3];
rz(-1.8281507) q[3];
sx q[3];
rz(-1.9338231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81994098) q[0];
sx q[0];
rz(-2.2714748) q[0];
sx q[0];
rz(-2.8735549) q[0];
rz(-0.23621121) q[1];
sx q[1];
rz(-1.7832489) q[1];
sx q[1];
rz(-1.4150298) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38430308) q[0];
sx q[0];
rz(-1.6727323) q[0];
sx q[0];
rz(-1.524964) q[0];
rz(-0.42534624) q[2];
sx q[2];
rz(-0.39452628) q[2];
sx q[2];
rz(0.91532842) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9795003) q[1];
sx q[1];
rz(-2.0071173) q[1];
sx q[1];
rz(1.8772582) q[1];
x q[2];
rz(2.799166) q[3];
sx q[3];
rz(-0.70243109) q[3];
sx q[3];
rz(0.46063761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1872824) q[2];
sx q[2];
rz(-2.8128746) q[2];
sx q[2];
rz(1.3817361) q[2];
rz(2.7749744) q[3];
sx q[3];
rz(-1.6691875) q[3];
sx q[3];
rz(0.097675145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5496552) q[0];
sx q[0];
rz(-2.9357935) q[0];
sx q[0];
rz(3.1254712) q[0];
rz(-1.4664117) q[1];
sx q[1];
rz(-1.5212395) q[1];
sx q[1];
rz(-0.88502562) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21986905) q[0];
sx q[0];
rz(-1.9737287) q[0];
sx q[0];
rz(-1.3833439) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22815223) q[2];
sx q[2];
rz(-1.9288262) q[2];
sx q[2];
rz(1.4756853) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.72674561) q[1];
sx q[1];
rz(-1.1903084) q[1];
sx q[1];
rz(0.49729113) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2161189) q[3];
sx q[3];
rz(-0.61788117) q[3];
sx q[3];
rz(0.91922744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.57899388) q[2];
sx q[2];
rz(-2.6851974) q[2];
sx q[2];
rz(-1.5835416) q[2];
rz(1.7700178) q[3];
sx q[3];
rz(-1.8669502) q[3];
sx q[3];
rz(0.19545999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1564932) q[0];
sx q[0];
rz(-1.3753966) q[0];
sx q[0];
rz(0.14955713) q[0];
rz(1.043148) q[1];
sx q[1];
rz(-0.97277343) q[1];
sx q[1];
rz(2.3057888) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0319894) q[0];
sx q[0];
rz(-1.2471203) q[0];
sx q[0];
rz(2.8780185) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.582242) q[2];
sx q[2];
rz(-0.69037333) q[2];
sx q[2];
rz(2.9857218) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0440825) q[1];
sx q[1];
rz(-2.9335576) q[1];
sx q[1];
rz(1.5745622) q[1];
x q[2];
rz(1.3204783) q[3];
sx q[3];
rz(-1.3930818) q[3];
sx q[3];
rz(1.5091997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0189556) q[2];
sx q[2];
rz(-1.476172) q[2];
sx q[2];
rz(2.7254851) q[2];
rz(-3.05919) q[3];
sx q[3];
rz(-0.96937886) q[3];
sx q[3];
rz(-2.940322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(3.0788954) q[0];
sx q[0];
rz(-1.0286999) q[0];
sx q[0];
rz(-1.0908352) q[0];
rz(2.003032) q[1];
sx q[1];
rz(-1.9416315) q[1];
sx q[1];
rz(0.60637766) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9370458) q[0];
sx q[0];
rz(-1.318745) q[0];
sx q[0];
rz(2.4840647) q[0];
rz(-0.80198432) q[2];
sx q[2];
rz(-1.9030266) q[2];
sx q[2];
rz(1.8301682) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7650827) q[1];
sx q[1];
rz(-2.174859) q[1];
sx q[1];
rz(2.2187697) q[1];
rz(0.71833937) q[3];
sx q[3];
rz(-1.6629617) q[3];
sx q[3];
rz(-0.80764333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.828317) q[2];
sx q[2];
rz(-0.54855359) q[2];
sx q[2];
rz(-0.11087785) q[2];
rz(-1.5781135) q[3];
sx q[3];
rz(-2.9302247) q[3];
sx q[3];
rz(1.3255239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.204708) q[0];
sx q[0];
rz(-0.81975833) q[0];
sx q[0];
rz(-2.4844266) q[0];
rz(-1.3232629) q[1];
sx q[1];
rz(-0.75138775) q[1];
sx q[1];
rz(-1.125186) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3599676) q[0];
sx q[0];
rz(-0.021718135) q[0];
sx q[0];
rz(2.3853073) q[0];
x q[1];
rz(-0.29203307) q[2];
sx q[2];
rz(-2.8148373) q[2];
sx q[2];
rz(2.9346443) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4954056) q[1];
sx q[1];
rz(-2.0584724) q[1];
sx q[1];
rz(1.133092) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7636307) q[3];
sx q[3];
rz(-0.71092793) q[3];
sx q[3];
rz(2.7148249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7323759) q[2];
sx q[2];
rz(-1.9595307) q[2];
sx q[2];
rz(2.1803975) q[2];
rz(3.1327278) q[3];
sx q[3];
rz(-0.88518393) q[3];
sx q[3];
rz(1.9302906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0601785) q[0];
sx q[0];
rz(-3.0280805) q[0];
sx q[0];
rz(0.42763448) q[0];
rz(1.112452) q[1];
sx q[1];
rz(-1.3149657) q[1];
sx q[1];
rz(2.1449259) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8450981) q[0];
sx q[0];
rz(-1.0834435) q[0];
sx q[0];
rz(-2.4250406) q[0];
rz(2.8585343) q[2];
sx q[2];
rz(-2.3583204) q[2];
sx q[2];
rz(-2.9976792) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.06293077) q[1];
sx q[1];
rz(-2.1306296) q[1];
sx q[1];
rz(-2.8275851) q[1];
rz(-pi) q[2];
rz(0.57076591) q[3];
sx q[3];
rz(-1.7460572) q[3];
sx q[3];
rz(1.4890763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4619649) q[2];
sx q[2];
rz(-2.6042016) q[2];
sx q[2];
rz(-2.4930387) q[2];
rz(0.31351659) q[3];
sx q[3];
rz(-2.2532012) q[3];
sx q[3];
rz(0.052481767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40470966) q[0];
sx q[0];
rz(-0.22060224) q[0];
sx q[0];
rz(-2.1746461) q[0];
rz(-0.95775682) q[1];
sx q[1];
rz(-2.0499332) q[1];
sx q[1];
rz(-1.0831833) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1304798) q[0];
sx q[0];
rz(-0.50043538) q[0];
sx q[0];
rz(2.6088663) q[0];
rz(-pi) q[1];
rz(-2.2317829) q[2];
sx q[2];
rz(-2.1860115) q[2];
sx q[2];
rz(1.3843975) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5965278) q[1];
sx q[1];
rz(-1.9633496) q[1];
sx q[1];
rz(2.3190581) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1414491) q[3];
sx q[3];
rz(-0.4793491) q[3];
sx q[3];
rz(1.0841313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0208685) q[2];
sx q[2];
rz(-2.3126297) q[2];
sx q[2];
rz(-1.4618358) q[2];
rz(0.36137897) q[3];
sx q[3];
rz(-1.8166108) q[3];
sx q[3];
rz(-0.98226515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7892889) q[0];
sx q[0];
rz(-2.4109349) q[0];
sx q[0];
rz(2.5590382) q[0];
rz(1.3652868) q[1];
sx q[1];
rz(-1.0368232) q[1];
sx q[1];
rz(0.22886151) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6270646) q[0];
sx q[0];
rz(-0.43958966) q[0];
sx q[0];
rz(2.702781) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9624356) q[2];
sx q[2];
rz(-1.8300042) q[2];
sx q[2];
rz(1.955206) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.072254213) q[1];
sx q[1];
rz(-1.7295341) q[1];
sx q[1];
rz(-0.065541288) q[1];
rz(-pi) q[2];
rz(1.1353605) q[3];
sx q[3];
rz(-1.2139911) q[3];
sx q[3];
rz(2.8716068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2710496) q[2];
sx q[2];
rz(-2.5898263) q[2];
sx q[2];
rz(-1.7331227) q[2];
rz(2.4208141) q[3];
sx q[3];
rz(-1.1495138) q[3];
sx q[3];
rz(-1.6225947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.085717399) q[0];
sx q[0];
rz(-1.7039104) q[0];
sx q[0];
rz(-1.3088551) q[0];
rz(-1.5720639) q[1];
sx q[1];
rz(-2.5253898) q[1];
sx q[1];
rz(0.22402221) q[1];
rz(-2.8926579) q[2];
sx q[2];
rz(-2.1215083) q[2];
sx q[2];
rz(-0.10100867) q[2];
rz(0.15275501) q[3];
sx q[3];
rz(-1.412231) q[3];
sx q[3];
rz(-2.9291736) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
