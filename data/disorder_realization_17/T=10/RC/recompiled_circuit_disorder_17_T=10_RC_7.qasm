OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.9085812) q[0];
sx q[0];
rz(4.3281875) q[0];
sx q[0];
rz(8.2789658) q[0];
rz(-2.1687578) q[1];
sx q[1];
rz(-1.6701148) q[1];
sx q[1];
rz(-0.29247984) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9777893) q[0];
sx q[0];
rz(-1.5153432) q[0];
sx q[0];
rz(-1.9652912) q[0];
rz(-pi) q[1];
rz(0.77941676) q[2];
sx q[2];
rz(-1.5343108) q[2];
sx q[2];
rz(1.2365637) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7293207) q[1];
sx q[1];
rz(-2.1643157) q[1];
sx q[1];
rz(0.055298474) q[1];
rz(-pi) q[2];
rz(2.0446159) q[3];
sx q[3];
rz(-1.5690656) q[3];
sx q[3];
rz(2.2497821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4108489) q[2];
sx q[2];
rz(-2.0480672) q[2];
sx q[2];
rz(2.6921819) q[2];
rz(0.64569008) q[3];
sx q[3];
rz(-2.460545) q[3];
sx q[3];
rz(1.2908363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9280424) q[0];
sx q[0];
rz(-2.2556861) q[0];
sx q[0];
rz(-0.74203062) q[0];
rz(1.6702601) q[1];
sx q[1];
rz(-0.69683087) q[1];
sx q[1];
rz(1.3630294) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5308967) q[0];
sx q[0];
rz(-0.96141059) q[0];
sx q[0];
rz(-1.0884398) q[0];
rz(-pi) q[1];
rz(1.3554553) q[2];
sx q[2];
rz(-2.0305579) q[2];
sx q[2];
rz(-0.41972566) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.14582536) q[1];
sx q[1];
rz(-1.2312504) q[1];
sx q[1];
rz(-1.550436) q[1];
x q[2];
rz(-0.95008534) q[3];
sx q[3];
rz(-0.32334298) q[3];
sx q[3];
rz(2.8107373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8403975) q[2];
sx q[2];
rz(-0.49673721) q[2];
sx q[2];
rz(1.8691241) q[2];
rz(-2.8043591) q[3];
sx q[3];
rz(-1.1703346) q[3];
sx q[3];
rz(0.01005323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067588016) q[0];
sx q[0];
rz(-2.8019866) q[0];
sx q[0];
rz(-2.8787676) q[0];
rz(-0.35573959) q[1];
sx q[1];
rz(-1.4590615) q[1];
sx q[1];
rz(1.9062818) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6806157) q[0];
sx q[0];
rz(-1.5810228) q[0];
sx q[0];
rz(3.1264722) q[0];
x q[1];
rz(1.4624865) q[2];
sx q[2];
rz(-2.8998067) q[2];
sx q[2];
rz(2.0568648) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9062412) q[1];
sx q[1];
rz(-1.4655349) q[1];
sx q[1];
rz(-0.32354849) q[1];
rz(-pi) q[2];
rz(-2.2658306) q[3];
sx q[3];
rz(-1.2998298) q[3];
sx q[3];
rz(-2.392644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0170903) q[2];
sx q[2];
rz(-1.3829145) q[2];
sx q[2];
rz(-0.18033218) q[2];
rz(2.5056433) q[3];
sx q[3];
rz(-1.162581) q[3];
sx q[3];
rz(-0.37160555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76698774) q[0];
sx q[0];
rz(-2.7174482) q[0];
sx q[0];
rz(-0.71907991) q[0];
rz(0.51302296) q[1];
sx q[1];
rz(-1.9388371) q[1];
sx q[1];
rz(-1.0346574) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35095222) q[0];
sx q[0];
rz(-1.6296903) q[0];
sx q[0];
rz(-1.590593) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.047495202) q[2];
sx q[2];
rz(-1.7058027) q[2];
sx q[2];
rz(-0.53802711) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.82995854) q[1];
sx q[1];
rz(-1.8533851) q[1];
sx q[1];
rz(-2.4294873) q[1];
x q[2];
rz(2.8095874) q[3];
sx q[3];
rz(-0.49626353) q[3];
sx q[3];
rz(-0.23976025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.93306142) q[2];
sx q[2];
rz(-2.8716817) q[2];
sx q[2];
rz(-2.9053524) q[2];
rz(1.158372) q[3];
sx q[3];
rz(-1.002243) q[3];
sx q[3];
rz(-2.2896144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1119969) q[0];
sx q[0];
rz(-2.1313666) q[0];
sx q[0];
rz(0.90233666) q[0];
rz(-0.92102712) q[1];
sx q[1];
rz(-2.5245456) q[1];
sx q[1];
rz(2.5193118) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8392004) q[0];
sx q[0];
rz(-1.7532187) q[0];
sx q[0];
rz(1.2829885) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.12139608) q[2];
sx q[2];
rz(-2.4897794) q[2];
sx q[2];
rz(0.35169841) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.2052853) q[1];
sx q[1];
rz(-2.949932) q[1];
sx q[1];
rz(-2.0602517) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4683687) q[3];
sx q[3];
rz(-0.58727467) q[3];
sx q[3];
rz(2.4622963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2174012) q[2];
sx q[2];
rz(-1.8175586) q[2];
sx q[2];
rz(-3.1203111) q[2];
rz(-1.4422669) q[3];
sx q[3];
rz(-0.71443021) q[3];
sx q[3];
rz(-2.2309979) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6570046) q[0];
sx q[0];
rz(-0.5963043) q[0];
sx q[0];
rz(0.044629991) q[0];
rz(-1.0021098) q[1];
sx q[1];
rz(-1.0005181) q[1];
sx q[1];
rz(0.034428509) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8196626) q[0];
sx q[0];
rz(-2.9557883) q[0];
sx q[0];
rz(1.4968605) q[0];
rz(-pi) q[1];
rz(2.0738515) q[2];
sx q[2];
rz(-1.6776553) q[2];
sx q[2];
rz(-1.6018794) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5163102) q[1];
sx q[1];
rz(-1.0370266) q[1];
sx q[1];
rz(0.37955243) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19272007) q[3];
sx q[3];
rz(-2.4270714) q[3];
sx q[3];
rz(2.2477828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.78952152) q[2];
sx q[2];
rz(-1.2723203) q[2];
sx q[2];
rz(-0.96674031) q[2];
rz(-0.012185193) q[3];
sx q[3];
rz(-2.2120078) q[3];
sx q[3];
rz(-0.10908443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0063342) q[0];
sx q[0];
rz(-2.874458) q[0];
sx q[0];
rz(0.99408856) q[0];
rz(-0.055123568) q[1];
sx q[1];
rz(-1.2722641) q[1];
sx q[1];
rz(0.73928839) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8871043) q[0];
sx q[0];
rz(-1.4758037) q[0];
sx q[0];
rz(-1.4332921) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0871068) q[2];
sx q[2];
rz(-1.0488044) q[2];
sx q[2];
rz(-0.99624485) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6608097) q[1];
sx q[1];
rz(-2.0936831) q[1];
sx q[1];
rz(-1.3516264) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2823295) q[3];
sx q[3];
rz(-0.97086421) q[3];
sx q[3];
rz(-2.695042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.58549515) q[2];
sx q[2];
rz(-1.7332417) q[2];
sx q[2];
rz(-2.5778256) q[2];
rz(-2.901315) q[3];
sx q[3];
rz(-2.6085745) q[3];
sx q[3];
rz(-1.3747922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4450842) q[0];
sx q[0];
rz(-0.17906469) q[0];
sx q[0];
rz(0.58404303) q[0];
rz(-0.42075992) q[1];
sx q[1];
rz(-2.1003484) q[1];
sx q[1];
rz(2.8631794) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7137372) q[0];
sx q[0];
rz(-2.5579778) q[0];
sx q[0];
rz(-1.7996656) q[0];
x q[1];
rz(-0.43012302) q[2];
sx q[2];
rz(-2.7134036) q[2];
sx q[2];
rz(2.6870514) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6555357) q[1];
sx q[1];
rz(-1.2206519) q[1];
sx q[1];
rz(-0.23305594) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1823468) q[3];
sx q[3];
rz(-2.5591345) q[3];
sx q[3];
rz(2.9782481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2723096) q[2];
sx q[2];
rz(-0.81988207) q[2];
sx q[2];
rz(2.753624) q[2];
rz(1.7913473) q[3];
sx q[3];
rz(-2.6781121) q[3];
sx q[3];
rz(-0.2750245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824654) q[0];
sx q[0];
rz(-0.21331856) q[0];
sx q[0];
rz(2.4097089) q[0];
rz(-0.16648079) q[1];
sx q[1];
rz(-0.44547588) q[1];
sx q[1];
rz(3.0678715) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026422231) q[0];
sx q[0];
rz(-1.1367043) q[0];
sx q[0];
rz(-0.11154453) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97754064) q[2];
sx q[2];
rz(-2.9974555) q[2];
sx q[2];
rz(-1.5365212) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8718865) q[1];
sx q[1];
rz(-1.5333813) q[1];
sx q[1];
rz(-2.5579631) q[1];
rz(-pi) q[2];
rz(-1.4599667) q[3];
sx q[3];
rz(-1.0004527) q[3];
sx q[3];
rz(3.0487206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59447294) q[2];
sx q[2];
rz(-2.9056845) q[2];
sx q[2];
rz(-2.7129042) q[2];
rz(1.3171014) q[3];
sx q[3];
rz(-1.7842112) q[3];
sx q[3];
rz(-1.2111506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4093032) q[0];
sx q[0];
rz(-1.6560873) q[0];
sx q[0];
rz(-1.0428585) q[0];
rz(-1.6304784) q[1];
sx q[1];
rz(-1.379456) q[1];
sx q[1];
rz(1.3483378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6207328) q[0];
sx q[0];
rz(-0.11611406) q[0];
sx q[0];
rz(-1.3470879) q[0];
rz(2.3073763) q[2];
sx q[2];
rz(-1.5110821) q[2];
sx q[2];
rz(0.63876736) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1765602) q[1];
sx q[1];
rz(-1.6415879) q[1];
sx q[1];
rz(-1.2464347) q[1];
rz(-pi) q[2];
rz(1.7616368) q[3];
sx q[3];
rz(-1.8075426) q[3];
sx q[3];
rz(-0.56425795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3412791) q[2];
sx q[2];
rz(-0.51395243) q[2];
sx q[2];
rz(3.0467765) q[2];
rz(-1.9684277) q[3];
sx q[3];
rz(-1.4011819) q[3];
sx q[3];
rz(2.9828984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.512758) q[0];
sx q[0];
rz(-0.47369581) q[0];
sx q[0];
rz(1.0439903) q[0];
rz(1.6013153) q[1];
sx q[1];
rz(-1.5834783) q[1];
sx q[1];
rz(-1.6316354) q[1];
rz(1.3576635) q[2];
sx q[2];
rz(-1.6518946) q[2];
sx q[2];
rz(1.2988731) q[2];
rz(2.1223162) q[3];
sx q[3];
rz(-2.301579) q[3];
sx q[3];
rz(-0.35001128) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
