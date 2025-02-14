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
rz(1.8546328) q[1];
sx q[1];
rz(-1.9998963) q[1];
sx q[1];
rz(0.72007522) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.274668) q[0];
sx q[0];
rz(-0.56827456) q[0];
sx q[0];
rz(2.3734762) q[0];
rz(2.949657) q[2];
sx q[2];
rz(-1.766381) q[2];
sx q[2];
rz(1.5337613) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.9125058) q[1];
sx q[1];
rz(-2.6709963) q[1];
sx q[1];
rz(2.5060593) q[1];
rz(-pi) q[2];
rz(2.9769865) q[3];
sx q[3];
rz(-0.96658488) q[3];
sx q[3];
rz(-0.34863499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2893452) q[2];
sx q[2];
rz(-2.2699558) q[2];
sx q[2];
rz(-0.80725011) q[2];
rz(1.6926258) q[3];
sx q[3];
rz(-2.2117386) q[3];
sx q[3];
rz(2.4248031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(-0.1194696) q[0];
sx q[0];
rz(-1.8417646) q[0];
sx q[0];
rz(1.3053373) q[0];
rz(2.3985825) q[1];
sx q[1];
rz(-0.65736714) q[1];
sx q[1];
rz(-1.2974799) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1851705) q[0];
sx q[0];
rz(-1.3885986) q[0];
sx q[0];
rz(-1.9745898) q[0];
rz(0.44746872) q[2];
sx q[2];
rz(-2.0966242) q[2];
sx q[2];
rz(-1.7350137) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4549437) q[1];
sx q[1];
rz(-1.2746069) q[1];
sx q[1];
rz(-2.185883) q[1];
rz(-pi) q[2];
rz(-2.8415806) q[3];
sx q[3];
rz(-0.37676806) q[3];
sx q[3];
rz(0.11480162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3798736) q[2];
sx q[2];
rz(-0.97731176) q[2];
sx q[2];
rz(1.0404111) q[2];
rz(0.31446332) q[3];
sx q[3];
rz(-1.9324666) q[3];
sx q[3];
rz(-1.8837121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8085025) q[0];
sx q[0];
rz(-2.1603778) q[0];
sx q[0];
rz(-1.8044949) q[0];
rz(2.5143276) q[1];
sx q[1];
rz(-1.7341055) q[1];
sx q[1];
rz(-1.5339392) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3298137) q[0];
sx q[0];
rz(-2.3372991) q[0];
sx q[0];
rz(0.95383184) q[0];
rz(2.5467403) q[2];
sx q[2];
rz(-0.50083465) q[2];
sx q[2];
rz(0.55933096) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.99974224) q[1];
sx q[1];
rz(-0.51020422) q[1];
sx q[1];
rz(-0.38533089) q[1];
rz(-pi) q[2];
rz(-0.19948761) q[3];
sx q[3];
rz(-1.1987682) q[3];
sx q[3];
rz(2.2147199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4072121) q[2];
sx q[2];
rz(-0.89443365) q[2];
sx q[2];
rz(-1.9341932) q[2];
rz(2.1842128) q[3];
sx q[3];
rz(-1.9060241) q[3];
sx q[3];
rz(-2.4522219) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44416881) q[0];
sx q[0];
rz(-0.95132315) q[0];
sx q[0];
rz(-0.091766894) q[0];
rz(-2.9755039) q[1];
sx q[1];
rz(-1.8548012) q[1];
sx q[1];
rz(-2.1002358) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0386376) q[0];
sx q[0];
rz(-0.79651662) q[0];
sx q[0];
rz(-0.060925555) q[0];
x q[1];
rz(1.9729593) q[2];
sx q[2];
rz(-2.3170174) q[2];
sx q[2];
rz(-0.74092016) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4554259) q[1];
sx q[1];
rz(-2.0435026) q[1];
sx q[1];
rz(0.03852877) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0001767) q[3];
sx q[3];
rz(-1.4849471) q[3];
sx q[3];
rz(0.72535634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.2287717) q[2];
sx q[2];
rz(-2.3232338) q[2];
sx q[2];
rz(-1.6010326) q[2];
rz(2.1238964) q[3];
sx q[3];
rz(-1.3295659) q[3];
sx q[3];
rz(1.8890007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8215264) q[0];
sx q[0];
rz(-1.3952661) q[0];
sx q[0];
rz(-0.23859247) q[0];
rz(-2.353031) q[1];
sx q[1];
rz(-0.32222727) q[1];
sx q[1];
rz(-0.90071789) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4960605) q[0];
sx q[0];
rz(-1.1405611) q[0];
sx q[0];
rz(-1.2853464) q[0];
rz(-1.1815669) q[2];
sx q[2];
rz(-1.6578362) q[2];
sx q[2];
rz(2.4554676) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.5821735) q[1];
sx q[1];
rz(-1.6985849) q[1];
sx q[1];
rz(-2.726597) q[1];
x q[2];
rz(2.3524893) q[3];
sx q[3];
rz(-0.89653095) q[3];
sx q[3];
rz(1.2885119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0939402) q[2];
sx q[2];
rz(-1.36146) q[2];
sx q[2];
rz(1.7729887) q[2];
rz(2.7641344) q[3];
sx q[3];
rz(-0.82337514) q[3];
sx q[3];
rz(1.9774168) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1222328) q[0];
sx q[0];
rz(-2.7074809) q[0];
sx q[0];
rz(2.9121616) q[0];
rz(2.5550628) q[1];
sx q[1];
rz(-2.6507381) q[1];
sx q[1];
rz(1.6612926) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23585715) q[0];
sx q[0];
rz(-0.93556306) q[0];
sx q[0];
rz(2.6978289) q[0];
rz(1.8312127) q[2];
sx q[2];
rz(-1.0170817) q[2];
sx q[2];
rz(1.2147724) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6203998) q[1];
sx q[1];
rz(-1.7971276) q[1];
sx q[1];
rz(1.625555) q[1];
rz(-pi) q[2];
rz(1.3020206) q[3];
sx q[3];
rz(-2.44938) q[3];
sx q[3];
rz(0.48066329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.71969879) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.557068) q[0];
sx q[0];
rz(-1.5140336) q[0];
sx q[0];
rz(2.6051482) q[0];
rz(-3.0570807) q[1];
sx q[1];
rz(-0.5785431) q[1];
sx q[1];
rz(-1.1892148) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3670334) q[0];
sx q[0];
rz(-1.4167455) q[0];
sx q[0];
rz(-0.94512248) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68185735) q[2];
sx q[2];
rz(-1.8769662) q[2];
sx q[2];
rz(-3.0628391) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7658733) q[1];
sx q[1];
rz(-0.60085591) q[1];
sx q[1];
rz(0.99581666) q[1];
rz(2.9042894) q[3];
sx q[3];
rz(-2.0198698) q[3];
sx q[3];
rz(-0.16332353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9540003) q[2];
sx q[2];
rz(-1.8686998) q[2];
sx q[2];
rz(-3.0583337) q[2];
rz(-1.8795053) q[3];
sx q[3];
rz(-0.83674651) q[3];
sx q[3];
rz(-2.4038103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5872203) q[0];
sx q[0];
rz(-1.6393336) q[0];
sx q[0];
rz(2.3578405) q[0];
rz(-1.6571677) q[1];
sx q[1];
rz(-2.617651) q[1];
sx q[1];
rz(-0.2934244) q[1];
rz(pi/2) q[2];
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
rz(-pi) q[1];
rz(1.8510455) q[2];
sx q[2];
rz(-2.4538605) q[2];
sx q[2];
rz(-2.6838145) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8415268) q[1];
sx q[1];
rz(-1.7151388) q[1];
sx q[1];
rz(-1.9886024) q[1];
x q[2];
rz(0.9154207) q[3];
sx q[3];
rz(-1.3815945) q[3];
sx q[3];
rz(2.9331741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36711127) q[2];
sx q[2];
rz(-1.4159055) q[2];
sx q[2];
rz(-1.6488546) q[2];
rz(-1.1501009) q[3];
sx q[3];
rz(-2.8425358) q[3];
sx q[3];
rz(2.2404631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8915326) q[0];
sx q[0];
rz(-1.0279259) q[0];
sx q[0];
rz(-3.133339) q[0];
rz(-0.0025657733) q[1];
sx q[1];
rz(-2.6305514) q[1];
sx q[1];
rz(2.0770238) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9998523) q[0];
sx q[0];
rz(-2.2069262) q[0];
sx q[0];
rz(-2.0881235) q[0];
rz(-pi) q[1];
rz(-3.1334115) q[2];
sx q[2];
rz(-1.3213833) q[2];
sx q[2];
rz(1.8298263) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9843205) q[1];
sx q[1];
rz(-2.0120579) q[1];
sx q[1];
rz(2.1391137) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58925924) q[3];
sx q[3];
rz(-1.4472826) q[3];
sx q[3];
rz(-0.46621399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0801487) q[2];
sx q[2];
rz(-1.6950357) q[2];
sx q[2];
rz(-2.9385369) q[2];
rz(-1.9301346) q[3];
sx q[3];
rz(-2.0720427) q[3];
sx q[3];
rz(3.0058461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9515297) q[0];
sx q[0];
rz(-2.0385346) q[0];
sx q[0];
rz(2.4679825) q[0];
rz(-2.8296962) q[1];
sx q[1];
rz(-1.2047647) q[1];
sx q[1];
rz(1.9685251) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86248428) q[0];
sx q[0];
rz(-2.0983834) q[0];
sx q[0];
rz(2.8010445) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3854957) q[2];
sx q[2];
rz(-1.7420766) q[2];
sx q[2];
rz(1.1403699) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0153326) q[1];
sx q[1];
rz(-1.646572) q[1];
sx q[1];
rz(-1.1557691) q[1];
rz(-pi) q[2];
rz(-2.2820149) q[3];
sx q[3];
rz(-1.2378581) q[3];
sx q[3];
rz(2.9200819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6761026) q[2];
sx q[2];
rz(-2.5238621) q[2];
sx q[2];
rz(0.22107302) q[2];
rz(2.5681833) q[3];
sx q[3];
rz(-2.2677877) q[3];
sx q[3];
rz(-2.0360086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7318646) q[0];
sx q[0];
rz(-1.8479713) q[0];
sx q[0];
rz(0.027298409) q[0];
rz(-1.9505386) q[1];
sx q[1];
rz(-1.7212894) q[1];
sx q[1];
rz(-1.7421834) q[1];
rz(1.2192192) q[2];
sx q[2];
rz(-0.58813358) q[2];
sx q[2];
rz(1.0434601) q[2];
rz(-0.97597573) q[3];
sx q[3];
rz(-1.2479758) q[3];
sx q[3];
rz(1.6231619) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
