OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9159311) q[0];
sx q[0];
rz(-0.8684648) q[0];
sx q[0];
rz(2.948569) q[0];
rz(-1.9999737) q[1];
sx q[1];
rz(3.5715754) q[1];
sx q[1];
rz(6.9663098) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8773168) q[0];
sx q[0];
rz(-1.5125934) q[0];
sx q[0];
rz(0.067406128) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97857742) q[2];
sx q[2];
rz(-0.41523146) q[2];
sx q[2];
rz(0.65939553) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4343623) q[1];
sx q[1];
rz(-1.1303567) q[1];
sx q[1];
rz(-0.67727725) q[1];
rz(-pi) q[2];
rz(0.37791032) q[3];
sx q[3];
rz(-1.0381191) q[3];
sx q[3];
rz(2.7045254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1258939) q[2];
sx q[2];
rz(-1.761972) q[2];
sx q[2];
rz(-2.0430298) q[2];
rz(-1.0788318) q[3];
sx q[3];
rz(-0.94509411) q[3];
sx q[3];
rz(-1.1014972) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1612448) q[0];
sx q[0];
rz(-0.315936) q[0];
sx q[0];
rz(0.20794491) q[0];
rz(-0.57693276) q[1];
sx q[1];
rz(-2.2536342) q[1];
sx q[1];
rz(1.6764486) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0315096) q[0];
sx q[0];
rz(-2.2950036) q[0];
sx q[0];
rz(1.5741041) q[0];
x q[1];
rz(-2.0589774) q[2];
sx q[2];
rz(-1.7980051) q[2];
sx q[2];
rz(-2.8607334) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0456902) q[1];
sx q[1];
rz(-1.1046788) q[1];
sx q[1];
rz(-0.63506605) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7705455) q[3];
sx q[3];
rz(-0.78073946) q[3];
sx q[3];
rz(2.3649529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0186105) q[2];
sx q[2];
rz(-0.98615042) q[2];
sx q[2];
rz(0.1097651) q[2];
rz(2.5189853) q[3];
sx q[3];
rz(-2.770335) q[3];
sx q[3];
rz(1.6842779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.1784172) q[0];
sx q[0];
rz(-0.85997471) q[0];
sx q[0];
rz(0.51613581) q[0];
rz(0.57488817) q[1];
sx q[1];
rz(-2.2153885) q[1];
sx q[1];
rz(-2.3410472) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.752906) q[0];
sx q[0];
rz(-1.9770925) q[0];
sx q[0];
rz(-2.9857062) q[0];
x q[1];
rz(2.4750701) q[2];
sx q[2];
rz(-2.1608673) q[2];
sx q[2];
rz(0.18596622) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.70558207) q[1];
sx q[1];
rz(-1.5758621) q[1];
sx q[1];
rz(-0.41298053) q[1];
x q[2];
rz(-2.27182) q[3];
sx q[3];
rz(-1.8043892) q[3];
sx q[3];
rz(1.9529395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4642554) q[2];
sx q[2];
rz(-0.3192454) q[2];
sx q[2];
rz(-1.3595954) q[2];
rz(0.96757403) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(-1.7165855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99252218) q[0];
sx q[0];
rz(-1.2620121) q[0];
sx q[0];
rz(2.676679) q[0];
rz(0.34856302) q[1];
sx q[1];
rz(-2.8788853) q[1];
sx q[1];
rz(-1.0850614) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5269055) q[0];
sx q[0];
rz(-2.5589716) q[0];
sx q[0];
rz(2.7132062) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8934946) q[2];
sx q[2];
rz(-1.2144226) q[2];
sx q[2];
rz(-0.46596371) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3241868) q[1];
sx q[1];
rz(-1.7122014) q[1];
sx q[1];
rz(-1.6721339) q[1];
rz(2.5892341) q[3];
sx q[3];
rz(-0.68867749) q[3];
sx q[3];
rz(-0.83113447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.00502914) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(0.68112779) q[2];
rz(-2.629771) q[3];
sx q[3];
rz(-2.8183283) q[3];
sx q[3];
rz(0.26369035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8624449) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(-1.408668) q[0];
rz(0.43235835) q[1];
sx q[1];
rz(-0.84195781) q[1];
sx q[1];
rz(0.98168215) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9156993) q[0];
sx q[0];
rz(-0.75076538) q[0];
sx q[0];
rz(2.4322926) q[0];
rz(0.9972516) q[2];
sx q[2];
rz(-2.855636) q[2];
sx q[2];
rz(-0.85948247) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9202068) q[1];
sx q[1];
rz(-1.1077987) q[1];
sx q[1];
rz(-1.7358001) q[1];
x q[2];
rz(-2.1682407) q[3];
sx q[3];
rz(-1.5095599) q[3];
sx q[3];
rz(0.51613584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0354054) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(0.48941082) q[2];
rz(-1.0148467) q[3];
sx q[3];
rz(-0.5265407) q[3];
sx q[3];
rz(1.813252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6569825) q[0];
sx q[0];
rz(-2.9841612) q[0];
sx q[0];
rz(0.37242517) q[0];
rz(1.8107481) q[1];
sx q[1];
rz(-2.1061888) q[1];
sx q[1];
rz(2.9763124) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3806992) q[0];
sx q[0];
rz(-1.5415511) q[0];
sx q[0];
rz(-3.0464892) q[0];
rz(-pi) q[1];
x q[1];
rz(0.68768244) q[2];
sx q[2];
rz(-1.5167987) q[2];
sx q[2];
rz(0.24535594) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2051516) q[1];
sx q[1];
rz(-1.1679822) q[1];
sx q[1];
rz(2.5564297) q[1];
rz(0.23541707) q[3];
sx q[3];
rz(-2.0093007) q[3];
sx q[3];
rz(-0.32270839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.91339397) q[2];
sx q[2];
rz(-1.0281576) q[2];
sx q[2];
rz(-1.6983263) q[2];
rz(-2.7741487) q[3];
sx q[3];
rz(-1.8668709) q[3];
sx q[3];
rz(0.2120367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3570324) q[0];
sx q[0];
rz(-2.050188) q[0];
sx q[0];
rz(-0.34926397) q[0];
rz(2.3941984) q[1];
sx q[1];
rz(-2.8458197) q[1];
sx q[1];
rz(-2.4051037) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38024494) q[0];
sx q[0];
rz(-1.5878829) q[0];
sx q[0];
rz(1.6197617) q[0];
rz(0.24918208) q[2];
sx q[2];
rz(-1.2417925) q[2];
sx q[2];
rz(2.5788384) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6474825) q[1];
sx q[1];
rz(-0.78299114) q[1];
sx q[1];
rz(-2.8857735) q[1];
rz(-pi) q[2];
rz(-0.50076671) q[3];
sx q[3];
rz(-1.9541249) q[3];
sx q[3];
rz(-2.4367743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13654576) q[2];
sx q[2];
rz(-1.9133291) q[2];
sx q[2];
rz(-2.6100256) q[2];
rz(0.68938869) q[3];
sx q[3];
rz(-1.6770984) q[3];
sx q[3];
rz(-1.2954856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.078995973) q[0];
sx q[0];
rz(-1.9364708) q[0];
sx q[0];
rz(1.8435562) q[0];
rz(0.80728665) q[1];
sx q[1];
rz(-1.1785945) q[1];
sx q[1];
rz(2.2198026) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73496504) q[0];
sx q[0];
rz(-1.7056744) q[0];
sx q[0];
rz(-0.69358967) q[0];
x q[1];
rz(0.83306649) q[2];
sx q[2];
rz(-1.7852011) q[2];
sx q[2];
rz(0.094878541) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3075879) q[1];
sx q[1];
rz(-1.3655791) q[1];
sx q[1];
rz(0.25968857) q[1];
rz(0.17955762) q[3];
sx q[3];
rz(-2.2710544) q[3];
sx q[3];
rz(1.2683271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.90199295) q[2];
sx q[2];
rz(-2.5805876) q[2];
sx q[2];
rz(-1.9699338) q[2];
rz(-1.7840067) q[3];
sx q[3];
rz(-1.4566908) q[3];
sx q[3];
rz(-1.4484891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25306025) q[0];
sx q[0];
rz(-1.1760412) q[0];
sx q[0];
rz(-1.4755479) q[0];
rz(-1.5400003) q[1];
sx q[1];
rz(-1.6801291) q[1];
sx q[1];
rz(0.67970651) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4363791) q[0];
sx q[0];
rz(-2.0853015) q[0];
sx q[0];
rz(1.3374469) q[0];
rz(0.18025132) q[2];
sx q[2];
rz(-0.88353523) q[2];
sx q[2];
rz(0.3790516) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3956086) q[1];
sx q[1];
rz(-1.9095699) q[1];
sx q[1];
rz(0.70396522) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1867832) q[3];
sx q[3];
rz(-0.98815742) q[3];
sx q[3];
rz(0.59347502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0150962) q[2];
sx q[2];
rz(-1.6588147) q[2];
sx q[2];
rz(-0.91040197) q[2];
rz(-0.67534584) q[3];
sx q[3];
rz(-0.93674913) q[3];
sx q[3];
rz(-2.2909686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05242059) q[0];
sx q[0];
rz(-1.2344673) q[0];
sx q[0];
rz(-0.7146548) q[0];
rz(0.71406281) q[1];
sx q[1];
rz(-0.95497447) q[1];
sx q[1];
rz(1.9649327) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.048621) q[0];
sx q[0];
rz(-1.7341359) q[0];
sx q[0];
rz(-0.049157814) q[0];
x q[1];
rz(-0.54729692) q[2];
sx q[2];
rz(-1.0101057) q[2];
sx q[2];
rz(2.1292357) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6224222) q[1];
sx q[1];
rz(-1.6455489) q[1];
sx q[1];
rz(-2.1291817) q[1];
rz(-pi) q[2];
x q[2];
rz(1.482974) q[3];
sx q[3];
rz(-1.0613958) q[3];
sx q[3];
rz(-2.1109964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.24361336) q[2];
sx q[2];
rz(-1.4695797) q[2];
sx q[2];
rz(2.0142377) q[2];
rz(-0.35774287) q[3];
sx q[3];
rz(-0.8711516) q[3];
sx q[3];
rz(-2.0991142) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42416278) q[0];
sx q[0];
rz(-1.2107727) q[0];
sx q[0];
rz(-2.6834224) q[0];
rz(-0.39623109) q[1];
sx q[1];
rz(-3.1165262) q[1];
sx q[1];
rz(-2.810626) q[1];
rz(-2.4049315) q[2];
sx q[2];
rz(-1.7792637) q[2];
sx q[2];
rz(-1.7532495) q[2];
rz(-0.022396537) q[3];
sx q[3];
rz(-2.7907484) q[3];
sx q[3];
rz(-0.69805935) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];