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
rz(0.32938862) q[0];
sx q[0];
rz(-2.5100799) q[0];
sx q[0];
rz(0.00087498571) q[0];
rz(-0.64088351) q[1];
sx q[1];
rz(5.2823245) q[1];
sx q[1];
rz(9.7699788) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.810173) q[0];
sx q[0];
rz(-1.6616482) q[0];
sx q[0];
rz(3.1167555) q[0];
x q[1];
rz(2.7617747) q[2];
sx q[2];
rz(-1.4336515) q[2];
sx q[2];
rz(-1.4963381) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6908603) q[1];
sx q[1];
rz(-1.1530515) q[1];
sx q[1];
rz(-0.91207204) q[1];
x q[2];
rz(-1.4963989) q[3];
sx q[3];
rz(-1.5870023) q[3];
sx q[3];
rz(2.5287573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.22902809) q[2];
sx q[2];
rz(-1.3338858) q[2];
sx q[2];
rz(-2.933617) q[2];
rz(2.8424756) q[3];
sx q[3];
rz(-0.57204539) q[3];
sx q[3];
rz(1.1408898) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8107373) q[0];
sx q[0];
rz(-2.9006697) q[0];
sx q[0];
rz(0.99288565) q[0];
rz(-1.317124) q[1];
sx q[1];
rz(-2.8254852) q[1];
sx q[1];
rz(-0.77450007) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68607722) q[0];
sx q[0];
rz(-2.9630399) q[0];
sx q[0];
rz(-1.5088085) q[0];
rz(-pi) q[1];
x q[1];
rz(2.108091) q[2];
sx q[2];
rz(-0.8391434) q[2];
sx q[2];
rz(3.0947859) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44821139) q[1];
sx q[1];
rz(-0.44084545) q[1];
sx q[1];
rz(0.54852672) q[1];
rz(1.9949739) q[3];
sx q[3];
rz(-0.50809233) q[3];
sx q[3];
rz(1.0095694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2449067) q[2];
sx q[2];
rz(-1.8550355) q[2];
sx q[2];
rz(-1.0150821) q[2];
rz(-0.24909881) q[3];
sx q[3];
rz(-2.2738012) q[3];
sx q[3];
rz(0.36756137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2772813) q[0];
sx q[0];
rz(-0.69121498) q[0];
sx q[0];
rz(0.15750289) q[0];
rz(1.0109673) q[1];
sx q[1];
rz(-0.33282655) q[1];
sx q[1];
rz(-0.13883042) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1942753) q[0];
sx q[0];
rz(-1.2091769) q[0];
sx q[0];
rz(0.34644923) q[0];
rz(1.9857668) q[2];
sx q[2];
rz(-0.89902821) q[2];
sx q[2];
rz(2.6443554) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5032387) q[1];
sx q[1];
rz(-1.279083) q[1];
sx q[1];
rz(-0.20511638) q[1];
x q[2];
rz(-2.8582358) q[3];
sx q[3];
rz(-1.4025926) q[3];
sx q[3];
rz(0.12780549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6802754) q[2];
sx q[2];
rz(-0.91247827) q[2];
sx q[2];
rz(0.3581363) q[2];
rz(-2.4628468) q[3];
sx q[3];
rz(-2.1923784) q[3];
sx q[3];
rz(2.1024735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0963652) q[0];
sx q[0];
rz(-0.97310936) q[0];
sx q[0];
rz(-0.3048234) q[0];
rz(2.7335956) q[1];
sx q[1];
rz(-1.7034737) q[1];
sx q[1];
rz(0.92794424) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.613425) q[0];
sx q[0];
rz(-1.4425462) q[0];
sx q[0];
rz(-1.7665461) q[0];
x q[1];
rz(-2.329823) q[2];
sx q[2];
rz(-0.78321811) q[2];
sx q[2];
rz(0.70657544) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.4907704) q[1];
sx q[1];
rz(-0.70413744) q[1];
sx q[1];
rz(0.23475523) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0545066) q[3];
sx q[3];
rz(-2.1355503) q[3];
sx q[3];
rz(1.3206583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9503595) q[2];
sx q[2];
rz(-2.5706036) q[2];
sx q[2];
rz(2.9534269) q[2];
rz(-1.865271) q[3];
sx q[3];
rz(-1.3151582) q[3];
sx q[3];
rz(1.4307384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4417878) q[0];
sx q[0];
rz(-0.34407523) q[0];
sx q[0];
rz(-0.019388327) q[0];
rz(-2.0806606) q[1];
sx q[1];
rz(-1.414199) q[1];
sx q[1];
rz(1.2394989) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65514046) q[0];
sx q[0];
rz(-1.0024655) q[0];
sx q[0];
rz(2.8080432) q[0];
rz(2.2197228) q[2];
sx q[2];
rz(-2.1989294) q[2];
sx q[2];
rz(-2.3482196) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.46249786) q[1];
sx q[1];
rz(-2.0789008) q[1];
sx q[1];
rz(-1.8930045) q[1];
x q[2];
rz(-2.0692213) q[3];
sx q[3];
rz(-1.2066505) q[3];
sx q[3];
rz(0.92723234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0197319) q[2];
sx q[2];
rz(-1.8283565) q[2];
sx q[2];
rz(1.3266374) q[2];
rz(-2.895368) q[3];
sx q[3];
rz(-1.1028057) q[3];
sx q[3];
rz(-2.2897913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5354079) q[0];
sx q[0];
rz(-0.98537213) q[0];
sx q[0];
rz(0.73915172) q[0];
rz(2.7958561) q[1];
sx q[1];
rz(-1.812499) q[1];
sx q[1];
rz(-0.87127042) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6046048) q[0];
sx q[0];
rz(-2.3571627) q[0];
sx q[0];
rz(-3.0434199) q[0];
x q[1];
rz(0.99314697) q[2];
sx q[2];
rz(-2.2340074) q[2];
sx q[2];
rz(-3.1301067) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.024684357) q[1];
sx q[1];
rz(-1.4977807) q[1];
sx q[1];
rz(-1.9597998) q[1];
rz(-0.63251791) q[3];
sx q[3];
rz(-0.7154724) q[3];
sx q[3];
rz(0.6882489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.31200108) q[2];
sx q[2];
rz(-2.5025554) q[2];
sx q[2];
rz(-0.87257067) q[2];
rz(1.5729337) q[3];
sx q[3];
rz(-0.60397732) q[3];
sx q[3];
rz(0.93305552) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0912112) q[0];
sx q[0];
rz(-0.25930431) q[0];
sx q[0];
rz(-0.76164371) q[0];
rz(-2.7161982) q[1];
sx q[1];
rz(-1.1401221) q[1];
sx q[1];
rz(0.92686191) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0935912) q[0];
sx q[0];
rz(-0.33080745) q[0];
sx q[0];
rz(0.76865102) q[0];
x q[1];
rz(0.92717391) q[2];
sx q[2];
rz(-0.96418751) q[2];
sx q[2];
rz(1.7546897) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.70434785) q[1];
sx q[1];
rz(-1.7926551) q[1];
sx q[1];
rz(1.1997486) q[1];
rz(-2.1939932) q[3];
sx q[3];
rz(-1.5316846) q[3];
sx q[3];
rz(0.87621237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1130134) q[2];
sx q[2];
rz(-2.4455652) q[2];
sx q[2];
rz(-0.87116233) q[2];
rz(0.29843676) q[3];
sx q[3];
rz(-1.8961743) q[3];
sx q[3];
rz(-1.3309853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7262064) q[0];
sx q[0];
rz(-0.50006777) q[0];
sx q[0];
rz(0.4183847) q[0];
rz(-0.2977953) q[1];
sx q[1];
rz(-1.1957542) q[1];
sx q[1];
rz(3.0072838) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90420049) q[0];
sx q[0];
rz(-0.89590329) q[0];
sx q[0];
rz(-2.2740721) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1654794) q[2];
sx q[2];
rz(-1.3994872) q[2];
sx q[2];
rz(-1.5754981) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1301535) q[1];
sx q[1];
rz(-1.7053889) q[1];
sx q[1];
rz(0.18421872) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28388309) q[3];
sx q[3];
rz(-1.1177269) q[3];
sx q[3];
rz(2.4456152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7316651) q[2];
sx q[2];
rz(-1.8525367) q[2];
sx q[2];
rz(0.64201391) q[2];
rz(2.0783453) q[3];
sx q[3];
rz(-0.65170538) q[3];
sx q[3];
rz(2.1003907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5409656) q[0];
sx q[0];
rz(-3.0525115) q[0];
sx q[0];
rz(-0.11216057) q[0];
rz(1.3345831) q[1];
sx q[1];
rz(-2.5490675) q[1];
sx q[1];
rz(-2.1122011) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6717259) q[0];
sx q[0];
rz(-0.70505667) q[0];
sx q[0];
rz(1.6204349) q[0];
rz(-pi) q[1];
x q[1];
rz(0.12597398) q[2];
sx q[2];
rz(-0.306189) q[2];
sx q[2];
rz(-2.9842215) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9255213) q[1];
sx q[1];
rz(-2.9072484) q[1];
sx q[1];
rz(-1.3804803) q[1];
rz(2.6736027) q[3];
sx q[3];
rz(-1.6525998) q[3];
sx q[3];
rz(-1.8379267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7182497) q[2];
sx q[2];
rz(-0.074904718) q[2];
sx q[2];
rz(0.038012803) q[2];
rz(2.5293317) q[3];
sx q[3];
rz(-0.93658787) q[3];
sx q[3];
rz(-0.14122252) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24899471) q[0];
sx q[0];
rz(-1.4709512) q[0];
sx q[0];
rz(-2.7227962) q[0];
rz(1.8839802) q[1];
sx q[1];
rz(-1.6323171) q[1];
sx q[1];
rz(-2.9990101) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2217956) q[0];
sx q[0];
rz(-0.16079535) q[0];
sx q[0];
rz(2.9662786) q[0];
rz(-pi) q[1];
rz(1.8480186) q[2];
sx q[2];
rz(-0.96053329) q[2];
sx q[2];
rz(-2.0137613) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0568472) q[1];
sx q[1];
rz(-2.1386302) q[1];
sx q[1];
rz(2.826087) q[1];
x q[2];
rz(0.94513388) q[3];
sx q[3];
rz(-0.54900733) q[3];
sx q[3];
rz(1.4956724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7964145) q[2];
sx q[2];
rz(-3.0609481) q[2];
sx q[2];
rz(-2.8734015) q[2];
rz(-1.4043407) q[3];
sx q[3];
rz(-1.0920478) q[3];
sx q[3];
rz(2.0973189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(1.1254697) q[0];
sx q[0];
rz(-0.37624993) q[0];
sx q[0];
rz(-2.7952623) q[0];
rz(3.012433) q[1];
sx q[1];
rz(-1.8864514) q[1];
sx q[1];
rz(1.4056978) q[1];
rz(-1.9228946) q[2];
sx q[2];
rz(-1.1209956) q[2];
sx q[2];
rz(1.6605177) q[2];
rz(0.066349647) q[3];
sx q[3];
rz(-1.2766311) q[3];
sx q[3];
rz(1.0286812) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
