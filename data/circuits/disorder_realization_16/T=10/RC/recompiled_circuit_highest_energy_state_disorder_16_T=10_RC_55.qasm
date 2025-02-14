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
rz(1.9749405) q[0];
sx q[0];
rz(-2.6986172) q[0];
sx q[0];
rz(-2.7583581) q[0];
rz(-1.2869599) q[1];
sx q[1];
rz(1.9998963) q[1];
sx q[1];
rz(13.286446) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0135277) q[0];
sx q[0];
rz(-1.1733455) q[0];
sx q[0];
rz(-1.1532408) q[0];
x q[1];
rz(-1.371648) q[2];
sx q[2];
rz(-1.3825644) q[2];
sx q[2];
rz(0.07478274) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2290869) q[1];
sx q[1];
rz(-2.6709963) q[1];
sx q[1];
rz(2.5060593) q[1];
x q[2];
rz(2.9769865) q[3];
sx q[3];
rz(-2.1750078) q[3];
sx q[3];
rz(0.34863499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2893452) q[2];
sx q[2];
rz(-0.87163681) q[2];
sx q[2];
rz(-0.80725011) q[2];
rz(-1.6926258) q[3];
sx q[3];
rz(-0.92985409) q[3];
sx q[3];
rz(-0.71678954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0221231) q[0];
sx q[0];
rz(-1.8417646) q[0];
sx q[0];
rz(-1.8362554) q[0];
rz(-0.7430101) q[1];
sx q[1];
rz(-0.65736714) q[1];
sx q[1];
rz(-1.2974799) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7867048) q[0];
sx q[0];
rz(-2.700665) q[0];
sx q[0];
rz(-1.132325) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99886151) q[2];
sx q[2];
rz(-1.9543658) q[2];
sx q[2];
rz(-3.0694196) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.686649) q[1];
sx q[1];
rz(-1.2746069) q[1];
sx q[1];
rz(-0.95570968) q[1];
x q[2];
rz(2.7801974) q[3];
sx q[3];
rz(-1.4618498) q[3];
sx q[3];
rz(1.7360842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.3798736) q[2];
sx q[2];
rz(-0.97731176) q[2];
sx q[2];
rz(1.0404111) q[2];
rz(2.8271293) q[3];
sx q[3];
rz(-1.209126) q[3];
sx q[3];
rz(1.2578806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.33309015) q[0];
sx q[0];
rz(-0.98121488) q[0];
sx q[0];
rz(-1.8044949) q[0];
rz(0.62726504) q[1];
sx q[1];
rz(-1.7341055) q[1];
sx q[1];
rz(1.5339392) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3298137) q[0];
sx q[0];
rz(-2.3372991) q[0];
sx q[0];
rz(-0.95383184) q[0];
x q[1];
rz(1.2731601) q[2];
sx q[2];
rz(-1.9797851) q[2];
sx q[2];
rz(-1.9252418) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9107289) q[1];
sx q[1];
rz(-1.3861935) q[1];
sx q[1];
rz(-0.47841013) q[1];
rz(-pi) q[2];
rz(0.19948761) q[3];
sx q[3];
rz(-1.1987682) q[3];
sx q[3];
rz(-2.2147199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.73438054) q[2];
sx q[2];
rz(-2.247159) q[2];
sx q[2];
rz(-1.9341932) q[2];
rz(-2.1842128) q[3];
sx q[3];
rz(-1.2355685) q[3];
sx q[3];
rz(0.68937075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(0.44416881) q[0];
sx q[0];
rz(-2.1902695) q[0];
sx q[0];
rz(-0.091766894) q[0];
rz(2.9755039) q[1];
sx q[1];
rz(-1.2867915) q[1];
sx q[1];
rz(-2.1002358) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4252151) q[0];
sx q[0];
rz(-1.6143403) q[0];
sx q[0];
rz(-2.3460044) q[0];
rz(-pi) q[1];
rz(1.1686334) q[2];
sx q[2];
rz(-0.82457525) q[2];
sx q[2];
rz(-0.74092016) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.68616679) q[1];
sx q[1];
rz(-1.09809) q[1];
sx q[1];
rz(0.03852877) q[1];
rz(0.54817537) q[3];
sx q[3];
rz(-0.16528567) q[3];
sx q[3];
rz(-0.30334869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.912821) q[2];
sx q[2];
rz(-0.81835881) q[2];
sx q[2];
rz(1.5405601) q[2];
rz(2.1238964) q[3];
sx q[3];
rz(-1.3295659) q[3];
sx q[3];
rz(1.8890007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8215264) q[0];
sx q[0];
rz(-1.3952661) q[0];
sx q[0];
rz(0.23859247) q[0];
rz(2.353031) q[1];
sx q[1];
rz(-0.32222727) q[1];
sx q[1];
rz(0.90071789) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.032388587) q[0];
sx q[0];
rz(-2.6302543) q[0];
sx q[0];
rz(-2.5912259) q[0];
rz(-pi) q[1];
rz(1.1815669) q[2];
sx q[2];
rz(-1.6578362) q[2];
sx q[2];
rz(0.6861251) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5594192) q[1];
sx q[1];
rz(-1.6985849) q[1];
sx q[1];
rz(0.41499563) q[1];
rz(0.78910335) q[3];
sx q[3];
rz(-0.89653095) q[3];
sx q[3];
rz(1.8530807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0939402) q[2];
sx q[2];
rz(-1.36146) q[2];
sx q[2];
rz(-1.7729887) q[2];
rz(-2.7641344) q[3];
sx q[3];
rz(-0.82337514) q[3];
sx q[3];
rz(-1.9774168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1222328) q[0];
sx q[0];
rz(-2.7074809) q[0];
sx q[0];
rz(-2.9121616) q[0];
rz(-2.5550628) q[1];
sx q[1];
rz(-0.49085453) q[1];
sx q[1];
rz(1.6612926) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.081588) q[0];
sx q[0];
rz(-1.2179273) q[0];
sx q[0];
rz(-2.2553483) q[0];
rz(-pi) q[1];
rz(0.39463674) q[2];
sx q[2];
rz(-0.60606128) q[2];
sx q[2];
rz(2.3958423) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0619032) q[1];
sx q[1];
rz(-1.5174355) q[1];
sx q[1];
rz(-0.22665955) q[1];
x q[2];
rz(1.3020206) q[3];
sx q[3];
rz(-2.44938) q[3];
sx q[3];
rz(-2.6609294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4218939) q[2];
sx q[2];
rz(-0.59013683) q[2];
sx q[2];
rz(1.3479007) q[2];
rz(-2.9501996) q[3];
sx q[3];
rz(-2.8189711) q[3];
sx q[3];
rz(1.9026683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58452463) q[0];
sx q[0];
rz(-1.6275591) q[0];
sx q[0];
rz(2.6051482) q[0];
rz(-3.0570807) q[1];
sx q[1];
rz(-2.5630496) q[1];
sx q[1];
rz(1.1892148) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7285648) q[0];
sx q[0];
rz(-2.4997093) q[0];
sx q[0];
rz(1.8299787) q[0];
rz(-pi) q[1];
rz(-0.68185735) q[2];
sx q[2];
rz(-1.2646265) q[2];
sx q[2];
rz(-3.0628391) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8457905) q[1];
sx q[1];
rz(-1.2582878) q[1];
sx q[1];
rz(1.048823) q[1];
rz(-pi) q[2];
rz(2.0246451) q[3];
sx q[3];
rz(-2.6374808) q[3];
sx q[3];
rz(-0.6716118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9540003) q[2];
sx q[2];
rz(-1.2728929) q[2];
sx q[2];
rz(3.0583337) q[2];
rz(1.2620874) q[3];
sx q[3];
rz(-2.3048461) q[3];
sx q[3];
rz(-0.73778233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5543723) q[0];
sx q[0];
rz(-1.502259) q[0];
sx q[0];
rz(-0.7837522) q[0];
rz(-1.6571677) q[1];
sx q[1];
rz(-2.617651) q[1];
sx q[1];
rz(-0.2934244) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83111887) q[0];
sx q[0];
rz(-1.9378512) q[0];
sx q[0];
rz(1.0793174) q[0];
rz(1.8510455) q[2];
sx q[2];
rz(-0.68773213) q[2];
sx q[2];
rz(2.6838145) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3345016) q[1];
sx q[1];
rz(-1.1576022) q[1];
sx q[1];
rz(0.15771247) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9154207) q[3];
sx q[3];
rz(-1.7599981) q[3];
sx q[3];
rz(-2.9331741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.36711127) q[2];
sx q[2];
rz(-1.4159055) q[2];
sx q[2];
rz(1.492738) q[2];
rz(1.1501009) q[3];
sx q[3];
rz(-0.29905683) q[3];
sx q[3];
rz(-0.90112954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25006008) q[0];
sx q[0];
rz(-2.1136668) q[0];
sx q[0];
rz(-3.133339) q[0];
rz(-0.0025657733) q[1];
sx q[1];
rz(-2.6305514) q[1];
sx q[1];
rz(-1.0645688) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24495793) q[0];
sx q[0];
rz(-1.9799398) q[0];
sx q[0];
rz(0.70434241) q[0];
rz(-pi) q[1];
rz(-1.6029036) q[2];
sx q[2];
rz(-0.24954441) q[2];
sx q[2];
rz(-1.3448993) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9656532) q[1];
sx q[1];
rz(-0.70427952) q[1];
sx q[1];
rz(-2.2910816) q[1];
rz(-2.5523334) q[3];
sx q[3];
rz(-1.6943101) q[3];
sx q[3];
rz(-2.6753787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0801487) q[2];
sx q[2];
rz(-1.4465569) q[2];
sx q[2];
rz(-0.20305571) q[2];
rz(-1.2114581) q[3];
sx q[3];
rz(-2.0720427) q[3];
sx q[3];
rz(-3.0058461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1900629) q[0];
sx q[0];
rz(-2.0385346) q[0];
sx q[0];
rz(0.67361012) q[0];
rz(2.8296962) q[1];
sx q[1];
rz(-1.936828) q[1];
sx q[1];
rz(1.9685251) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6098227) q[0];
sx q[0];
rz(-1.2780452) q[0];
sx q[0];
rz(-1.0171153) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3373702) q[2];
sx q[2];
rz(-2.3131822) q[2];
sx q[2];
rz(2.5517922) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7562779) q[1];
sx q[1];
rz(-0.42149252) q[1];
sx q[1];
rz(1.384686) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7133663) q[3];
sx q[3];
rz(-0.90598327) q[3];
sx q[3];
rz(1.5178102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46549002) q[2];
sx q[2];
rz(-0.61773053) q[2];
sx q[2];
rz(-0.22107302) q[2];
rz(-2.5681833) q[3];
sx q[3];
rz(-2.2677877) q[3];
sx q[3];
rz(-1.1055841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.7318646) q[0];
sx q[0];
rz(-1.2936214) q[0];
sx q[0];
rz(-3.1142942) q[0];
rz(-1.9505386) q[1];
sx q[1];
rz(-1.7212894) q[1];
sx q[1];
rz(-1.7421834) q[1];
rz(-1.0114317) q[2];
sx q[2];
rz(-1.7630429) q[2];
sx q[2];
rz(-0.2311308) q[2];
rz(2.1656169) q[3];
sx q[3];
rz(-1.2479758) q[3];
sx q[3];
rz(1.6231619) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
