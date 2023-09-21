OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4361753) q[0];
sx q[0];
rz(-0.56646148) q[0];
sx q[0];
rz(0.17106549) q[0];
rz(-1.6879727) q[1];
sx q[1];
rz(3.4847335) q[1];
sx q[1];
rz(7.6141678) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91141191) q[0];
sx q[0];
rz(-2.4898306) q[0];
sx q[0];
rz(-0.08641152) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51585977) q[2];
sx q[2];
rz(-1.5896279) q[2];
sx q[2];
rz(3.0992103) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6273856) q[1];
sx q[1];
rz(-1.5878907) q[1];
sx q[1];
rz(1.3196726) q[1];
x q[2];
rz(-0.29403789) q[3];
sx q[3];
rz(-0.66411823) q[3];
sx q[3];
rz(-2.3703863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.77582899) q[2];
sx q[2];
rz(-2.3684431) q[2];
sx q[2];
rz(1.3295056) q[2];
rz(-1.4154411) q[3];
sx q[3];
rz(-1.5664145) q[3];
sx q[3];
rz(-2.9392488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99130327) q[0];
sx q[0];
rz(-2.5973899) q[0];
sx q[0];
rz(-1.6888899) q[0];
rz(2.9406722) q[1];
sx q[1];
rz(-2.0566437) q[1];
sx q[1];
rz(2.8413049) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6712924) q[0];
sx q[0];
rz(-1.6912582) q[0];
sx q[0];
rz(0.60351535) q[0];
rz(-0.57057256) q[2];
sx q[2];
rz(-1.8405481) q[2];
sx q[2];
rz(2.3938993) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.86094942) q[1];
sx q[1];
rz(-1.7935392) q[1];
sx q[1];
rz(1.8722948) q[1];
x q[2];
rz(-2.2829451) q[3];
sx q[3];
rz(-2.894069) q[3];
sx q[3];
rz(1.991322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77256569) q[2];
sx q[2];
rz(-0.3521266) q[2];
sx q[2];
rz(2.1035813) q[2];
rz(2.299262) q[3];
sx q[3];
rz(-1.7269644) q[3];
sx q[3];
rz(-0.37030181) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5791941) q[0];
sx q[0];
rz(-2.6694522) q[0];
sx q[0];
rz(2.753479) q[0];
rz(3.0691222) q[1];
sx q[1];
rz(-1.4265172) q[1];
sx q[1];
rz(-0.31633502) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2280699) q[0];
sx q[0];
rz(-1.5689578) q[0];
sx q[0];
rz(-1.3257922) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5306273) q[2];
sx q[2];
rz(-2.4749711) q[2];
sx q[2];
rz(-1.9463469) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0132732) q[1];
sx q[1];
rz(-1.1228021) q[1];
sx q[1];
rz(-2.0112579) q[1];
rz(0.29052492) q[3];
sx q[3];
rz(-1.8968664) q[3];
sx q[3];
rz(2.3323004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1091653) q[2];
sx q[2];
rz(-0.18523231) q[2];
sx q[2];
rz(-0.16080984) q[2];
rz(-0.12604776) q[3];
sx q[3];
rz(-1.7536609) q[3];
sx q[3];
rz(-1.0236615) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2407103) q[0];
sx q[0];
rz(-0.70403376) q[0];
sx q[0];
rz(-0.8316935) q[0];
rz(1.3446993) q[1];
sx q[1];
rz(-0.99935499) q[1];
sx q[1];
rz(-1.5136738) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31878372) q[0];
sx q[0];
rz(-1.6167287) q[0];
sx q[0];
rz(2.4687693) q[0];
rz(-0.081110031) q[2];
sx q[2];
rz(-1.0081648) q[2];
sx q[2];
rz(-2.5155591) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0107683) q[1];
sx q[1];
rz(-2.0190034) q[1];
sx q[1];
rz(-3.1229449) q[1];
rz(-pi) q[2];
rz(1.1553331) q[3];
sx q[3];
rz(-2.2586125) q[3];
sx q[3];
rz(-2.8727939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4251129) q[2];
sx q[2];
rz(-2.1322865) q[2];
sx q[2];
rz(-1.224219) q[2];
rz(-2.7246357) q[3];
sx q[3];
rz(-2.0988393) q[3];
sx q[3];
rz(-0.52880374) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.058218) q[0];
sx q[0];
rz(-0.26979065) q[0];
sx q[0];
rz(1.8378687) q[0];
rz(2.6858792) q[1];
sx q[1];
rz(-2.8790751) q[1];
sx q[1];
rz(0.051503332) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18729678) q[0];
sx q[0];
rz(-1.4547537) q[0];
sx q[0];
rz(-0.10033484) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8565882) q[2];
sx q[2];
rz(-1.4225866) q[2];
sx q[2];
rz(-2.1697901) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3225973) q[1];
sx q[1];
rz(-1.7261191) q[1];
sx q[1];
rz(1.1197234) q[1];
rz(-1.6938985) q[3];
sx q[3];
rz(-0.8408635) q[3];
sx q[3];
rz(-0.52809944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.45433989) q[2];
sx q[2];
rz(-1.3663102) q[2];
sx q[2];
rz(-0.59801897) q[2];
rz(-2.5915742) q[3];
sx q[3];
rz(-2.9019182) q[3];
sx q[3];
rz(2.0050744) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0710058) q[0];
sx q[0];
rz(-3.0650009) q[0];
sx q[0];
rz(2.9396074) q[0];
rz(0.96549353) q[1];
sx q[1];
rz(-2.1172724) q[1];
sx q[1];
rz(3.0583256) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4250454) q[0];
sx q[0];
rz(-1.7444567) q[0];
sx q[0];
rz(2.9437149) q[0];
rz(-pi) q[1];
rz(-0.27481885) q[2];
sx q[2];
rz(-1.3126557) q[2];
sx q[2];
rz(-1.2061662) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.49230089) q[1];
sx q[1];
rz(-1.9363083) q[1];
sx q[1];
rz(3.027012) q[1];
rz(-pi) q[2];
rz(1.2721328) q[3];
sx q[3];
rz(-2.1080058) q[3];
sx q[3];
rz(2.1512973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.10963708) q[2];
sx q[2];
rz(-1.6836616) q[2];
sx q[2];
rz(1.7588245) q[2];
rz(-2.8921195) q[3];
sx q[3];
rz(-1.8981257) q[3];
sx q[3];
rz(-1.4613387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4270585) q[0];
sx q[0];
rz(-0.80474168) q[0];
sx q[0];
rz(0.89685857) q[0];
rz(-0.25009051) q[1];
sx q[1];
rz(-0.18083328) q[1];
sx q[1];
rz(1.1175964) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9761336) q[0];
sx q[0];
rz(-1.4384067) q[0];
sx q[0];
rz(-0.73029851) q[0];
rz(-pi) q[1];
rz(-2.2824077) q[2];
sx q[2];
rz(-1.4118328) q[2];
sx q[2];
rz(-2.8189916) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7270131) q[1];
sx q[1];
rz(-2.1257183) q[1];
sx q[1];
rz(-1.3347866) q[1];
rz(-pi) q[2];
rz(0.44495961) q[3];
sx q[3];
rz(-0.58044725) q[3];
sx q[3];
rz(-2.2821033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.08427944) q[2];
sx q[2];
rz(-1.7417615) q[2];
sx q[2];
rz(-2.9220667) q[2];
rz(-2.7548742) q[3];
sx q[3];
rz(-2.3733449) q[3];
sx q[3];
rz(0.99532551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052208386) q[0];
sx q[0];
rz(-1.5445671) q[0];
sx q[0];
rz(0.21729939) q[0];
rz(0.030933881) q[1];
sx q[1];
rz(-2.5035281) q[1];
sx q[1];
rz(-0.6589748) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12595183) q[0];
sx q[0];
rz(-0.77227393) q[0];
sx q[0];
rz(0.036235972) q[0];
rz(-0.89491567) q[2];
sx q[2];
rz(-2.1721828) q[2];
sx q[2];
rz(2.2611157) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.36775667) q[1];
sx q[1];
rz(-2.3564597) q[1];
sx q[1];
rz(-2.8326041) q[1];
x q[2];
rz(3.072776) q[3];
sx q[3];
rz(-2.2127082) q[3];
sx q[3];
rz(-1.9300234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7028246) q[2];
sx q[2];
rz(-1.7610901) q[2];
sx q[2];
rz(-0.32021114) q[2];
rz(-1.5252339) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(-1.2493791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3808688) q[0];
sx q[0];
rz(-0.18146935) q[0];
sx q[0];
rz(-0.6828126) q[0];
rz(2.071351) q[1];
sx q[1];
rz(-1.8990592) q[1];
sx q[1];
rz(-1.210093) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72737981) q[0];
sx q[0];
rz(-2.7276037) q[0];
sx q[0];
rz(-1.0340286) q[0];
rz(-pi) q[1];
rz(-2.8898583) q[2];
sx q[2];
rz(-2.1605957) q[2];
sx q[2];
rz(-2.7272607) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.71483892) q[1];
sx q[1];
rz(-2.4086082) q[1];
sx q[1];
rz(-2.0183802) q[1];
x q[2];
rz(2.5217767) q[3];
sx q[3];
rz(-1.2634988) q[3];
sx q[3];
rz(2.4710771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5658297) q[2];
sx q[2];
rz(-1.016022) q[2];
sx q[2];
rz(0.56662095) q[2];
rz(-2.2120655) q[3];
sx q[3];
rz(-0.98171392) q[3];
sx q[3];
rz(-1.367759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41314769) q[0];
sx q[0];
rz(-0.25756535) q[0];
sx q[0];
rz(-1.7096827) q[0];
rz(-0.52629772) q[1];
sx q[1];
rz(-0.52733517) q[1];
sx q[1];
rz(-2.3419103) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.66946) q[0];
sx q[0];
rz(-1.9181983) q[0];
sx q[0];
rz(3.0513289) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4893555) q[2];
sx q[2];
rz(-1.4105721) q[2];
sx q[2];
rz(-0.83934957) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.77196808) q[1];
sx q[1];
rz(-2.5110285) q[1];
sx q[1];
rz(0.97009138) q[1];
x q[2];
rz(-1.6465854) q[3];
sx q[3];
rz(-1.7160551) q[3];
sx q[3];
rz(2.4491572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8563103) q[2];
sx q[2];
rz(-2.6976863) q[2];
sx q[2];
rz(0.69520673) q[2];
rz(-0.55784145) q[3];
sx q[3];
rz(-1.7772243) q[3];
sx q[3];
rz(-0.69303304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0556864) q[0];
sx q[0];
rz(-1.9491371) q[0];
sx q[0];
rz(0.51167713) q[0];
rz(-1.862539) q[1];
sx q[1];
rz(-0.74395724) q[1];
sx q[1];
rz(-0.67768135) q[1];
rz(1.5699301) q[2];
sx q[2];
rz(-1.2348839) q[2];
sx q[2];
rz(0.15765794) q[2];
rz(-0.79808509) q[3];
sx q[3];
rz(-1.6860262) q[3];
sx q[3];
rz(-1.891914) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];