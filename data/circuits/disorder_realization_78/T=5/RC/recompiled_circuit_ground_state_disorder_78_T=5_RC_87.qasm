OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8371589) q[0];
sx q[0];
rz(4.9871939) q[0];
sx q[0];
rz(11.468588) q[0];
rz(-0.9552362) q[1];
sx q[1];
rz(-2.3999441) q[1];
sx q[1];
rz(-0.15129605) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3442952) q[0];
sx q[0];
rz(-1.5963285) q[0];
sx q[0];
rz(-1.4454557) q[0];
rz(-pi) q[1];
x q[1];
rz(2.120235) q[2];
sx q[2];
rz(-1.6701227) q[2];
sx q[2];
rz(1.8252357) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9469493) q[1];
sx q[1];
rz(-2.2990755) q[1];
sx q[1];
rz(0.76232736) q[1];
rz(-pi) q[2];
rz(1.8586848) q[3];
sx q[3];
rz(-2.3514053) q[3];
sx q[3];
rz(-1.1367281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1260881) q[2];
sx q[2];
rz(-1.8512923) q[2];
sx q[2];
rz(-2.8224831) q[2];
rz(1.1110405) q[3];
sx q[3];
rz(-2.5824661) q[3];
sx q[3];
rz(1.3148974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023271712) q[0];
sx q[0];
rz(-2.6225704) q[0];
sx q[0];
rz(0.58854377) q[0];
rz(-2.5449246) q[1];
sx q[1];
rz(-1.8110954) q[1];
sx q[1];
rz(-2.9002424) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4174735) q[0];
sx q[0];
rz(-0.80097526) q[0];
sx q[0];
rz(0.95306122) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0701837) q[2];
sx q[2];
rz(-1.6334849) q[2];
sx q[2];
rz(-0.78298616) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4171364) q[1];
sx q[1];
rz(-2.0057185) q[1];
sx q[1];
rz(0.10970727) q[1];
rz(-pi) q[2];
rz(-0.94542687) q[3];
sx q[3];
rz(-1.6540048) q[3];
sx q[3];
rz(2.0997467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1179463) q[2];
sx q[2];
rz(-1.4586552) q[2];
sx q[2];
rz(-0.092197593) q[2];
rz(0.89208952) q[3];
sx q[3];
rz(-2.2690319) q[3];
sx q[3];
rz(-2.0766855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43111619) q[0];
sx q[0];
rz(-1.5065864) q[0];
sx q[0];
rz(3.0431252) q[0];
rz(0.59421986) q[1];
sx q[1];
rz(-1.0585982) q[1];
sx q[1];
rz(2.1509511) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38483175) q[0];
sx q[0];
rz(-0.38479003) q[0];
sx q[0];
rz(-1.8033474) q[0];
x q[1];
rz(2.5712396) q[2];
sx q[2];
rz(-0.84707558) q[2];
sx q[2];
rz(-3.0511659) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.65730598) q[1];
sx q[1];
rz(-0.68731824) q[1];
sx q[1];
rz(1.1794075) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0001171) q[3];
sx q[3];
rz(-0.50715441) q[3];
sx q[3];
rz(2.0961268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.46382612) q[2];
sx q[2];
rz(-2.3550484) q[2];
sx q[2];
rz(2.3393935) q[2];
rz(-1.5944611) q[3];
sx q[3];
rz(-2.0764669) q[3];
sx q[3];
rz(-2.2772363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39740729) q[0];
sx q[0];
rz(-2.7825401) q[0];
sx q[0];
rz(1.3209976) q[0];
rz(1.6162704) q[1];
sx q[1];
rz(-2.1883712) q[1];
sx q[1];
rz(-0.1604518) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9994151) q[0];
sx q[0];
rz(-0.83926187) q[0];
sx q[0];
rz(2.8150303) q[0];
rz(-0.17244417) q[2];
sx q[2];
rz(-1.4566696) q[2];
sx q[2];
rz(0.44140377) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7509979) q[1];
sx q[1];
rz(-2.2143151) q[1];
sx q[1];
rz(-2.949763) q[1];
x q[2];
rz(0.64063756) q[3];
sx q[3];
rz(-0.63797073) q[3];
sx q[3];
rz(-0.20221329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.82115951) q[2];
sx q[2];
rz(-1.6511788) q[2];
sx q[2];
rz(2.5725345) q[2];
rz(-1.2218366) q[3];
sx q[3];
rz(-2.8025083) q[3];
sx q[3];
rz(-3.0088185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64567599) q[0];
sx q[0];
rz(-2.3585632) q[0];
sx q[0];
rz(-0.83086479) q[0];
rz(1.2105385) q[1];
sx q[1];
rz(-1.4930864) q[1];
sx q[1];
rz(2.1499706) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9308335) q[0];
sx q[0];
rz(-2.2607339) q[0];
sx q[0];
rz(0.4042952) q[0];
rz(-pi) q[1];
rz(0.28056552) q[2];
sx q[2];
rz(-0.90350752) q[2];
sx q[2];
rz(1.9633121) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2048671) q[1];
sx q[1];
rz(-0.74581742) q[1];
sx q[1];
rz(1.3406483) q[1];
rz(-pi) q[2];
rz(-2.1915221) q[3];
sx q[3];
rz(-2.3394428) q[3];
sx q[3];
rz(2.7190894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3811938) q[2];
sx q[2];
rz(-2.2187967) q[2];
sx q[2];
rz(2.7823616) q[2];
rz(0.71581101) q[3];
sx q[3];
rz(-2.3575213) q[3];
sx q[3];
rz(1.951096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0642218) q[0];
sx q[0];
rz(-1.2718028) q[0];
sx q[0];
rz(2.6203058) q[0];
rz(-0.84292665) q[1];
sx q[1];
rz(-1.9631674) q[1];
sx q[1];
rz(3.0311323) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.045521988) q[0];
sx q[0];
rz(-1.318249) q[0];
sx q[0];
rz(0.70742328) q[0];
rz(-1.6104389) q[2];
sx q[2];
rz(-2.4986598) q[2];
sx q[2];
rz(0.55243353) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1350491) q[1];
sx q[1];
rz(-1.1056756) q[1];
sx q[1];
rz(-2.9198482) q[1];
rz(1.0486097) q[3];
sx q[3];
rz(-1.3449826) q[3];
sx q[3];
rz(1.9584394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9342186) q[2];
sx q[2];
rz(-1.5997581) q[2];
sx q[2];
rz(0.99384394) q[2];
rz(-0.55073109) q[3];
sx q[3];
rz(-0.5144853) q[3];
sx q[3];
rz(1.5229092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9024502) q[0];
sx q[0];
rz(-2.7681594) q[0];
sx q[0];
rz(0.19110876) q[0];
rz(-0.36901078) q[1];
sx q[1];
rz(-1.7419107) q[1];
sx q[1];
rz(-0.63327995) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8275237) q[0];
sx q[0];
rz(-1.6681328) q[0];
sx q[0];
rz(0.23989664) q[0];
rz(-pi) q[1];
x q[1];
rz(0.52187829) q[2];
sx q[2];
rz(-0.79127914) q[2];
sx q[2];
rz(-2.996614) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.889588) q[1];
sx q[1];
rz(-1.377863) q[1];
sx q[1];
rz(-1.361752) q[1];
x q[2];
rz(-2.4426887) q[3];
sx q[3];
rz(-1.6972542) q[3];
sx q[3];
rz(-1.5918283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.17948991) q[2];
sx q[2];
rz(-1.7682163) q[2];
sx q[2];
rz(-2.7607259) q[2];
rz(1.7535836) q[3];
sx q[3];
rz(-2.1151147) q[3];
sx q[3];
rz(-1.4306205) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66985828) q[0];
sx q[0];
rz(-0.40100455) q[0];
sx q[0];
rz(1.5420472) q[0];
rz(-1.764864) q[1];
sx q[1];
rz(-1.7574666) q[1];
sx q[1];
rz(-0.9309887) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41808266) q[0];
sx q[0];
rz(-0.93659725) q[0];
sx q[0];
rz(-1.4618327) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3190218) q[2];
sx q[2];
rz(-1.1295756) q[2];
sx q[2];
rz(0.069442858) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.763995) q[1];
sx q[1];
rz(-0.74255172) q[1];
sx q[1];
rz(-1.5686036) q[1];
rz(-pi) q[2];
rz(-2.6205711) q[3];
sx q[3];
rz(-0.98257321) q[3];
sx q[3];
rz(-0.38834342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.57637438) q[2];
sx q[2];
rz(-0.56052506) q[2];
sx q[2];
rz(-3.0726688) q[2];
rz(-1.2636412) q[3];
sx q[3];
rz(-0.37216035) q[3];
sx q[3];
rz(2.4431958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7428335) q[0];
sx q[0];
rz(-0.95785207) q[0];
sx q[0];
rz(3.0897019) q[0];
rz(0.17008153) q[1];
sx q[1];
rz(-0.64337987) q[1];
sx q[1];
rz(0.35194078) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1387685) q[0];
sx q[0];
rz(-1.616298) q[0];
sx q[0];
rz(2.7753434) q[0];
x q[1];
rz(0.76100112) q[2];
sx q[2];
rz(-0.78520757) q[2];
sx q[2];
rz(-2.1987555) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4008949) q[1];
sx q[1];
rz(-2.4603421) q[1];
sx q[1];
rz(0.91886144) q[1];
rz(1.7991583) q[3];
sx q[3];
rz(-0.31410892) q[3];
sx q[3];
rz(2.4921162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.45741442) q[2];
sx q[2];
rz(-0.64932051) q[2];
sx q[2];
rz(-2.7049098) q[2];
rz(-0.39143482) q[3];
sx q[3];
rz(-1.6792363) q[3];
sx q[3];
rz(2.2552538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6201685) q[0];
sx q[0];
rz(-0.7069718) q[0];
sx q[0];
rz(-3.022505) q[0];
rz(1.8424312) q[1];
sx q[1];
rz(-1.3928587) q[1];
sx q[1];
rz(1.3577168) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8004476) q[0];
sx q[0];
rz(-0.94609944) q[0];
sx q[0];
rz(0.50826061) q[0];
rz(1.18612) q[2];
sx q[2];
rz(-2.3529077) q[2];
sx q[2];
rz(2.2269611) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.75666243) q[1];
sx q[1];
rz(-2.1637056) q[1];
sx q[1];
rz(1.97744) q[1];
rz(-2.1351027) q[3];
sx q[3];
rz(-2.3052633) q[3];
sx q[3];
rz(1.8759954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.208821) q[2];
sx q[2];
rz(-1.5522771) q[2];
sx q[2];
rz(-2.5271752) q[2];
rz(-1.6299853) q[3];
sx q[3];
rz(-2.4698518) q[3];
sx q[3];
rz(-0.083812788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83196249) q[0];
sx q[0];
rz(-0.93292581) q[0];
sx q[0];
rz(0.25136872) q[0];
rz(1.0328737) q[1];
sx q[1];
rz(-1.1943457) q[1];
sx q[1];
rz(1.7051382) q[1];
rz(-0.070214179) q[2];
sx q[2];
rz(-1.4660942) q[2];
sx q[2];
rz(1.6906665) q[2];
rz(0.80794215) q[3];
sx q[3];
rz(-2.1389037) q[3];
sx q[3];
rz(1.1566333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
