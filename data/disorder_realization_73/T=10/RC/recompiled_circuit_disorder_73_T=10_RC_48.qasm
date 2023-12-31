OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1459382) q[0];
sx q[0];
rz(3.6448195) q[0];
sx q[0];
rz(10.148944) q[0];
rz(0.63996285) q[1];
sx q[1];
rz(-0.53007403) q[1];
sx q[1];
rz(-0.78483265) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1544224) q[0];
sx q[0];
rz(-1.3599456) q[0];
sx q[0];
rz(0.2984557) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9169693) q[2];
sx q[2];
rz(-2.7135239) q[2];
sx q[2];
rz(-0.12878865) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6341056) q[1];
sx q[1];
rz(-1.7934985) q[1];
sx q[1];
rz(1.0313862) q[1];
rz(-pi) q[2];
rz(-1.4943487) q[3];
sx q[3];
rz(-2.7259698) q[3];
sx q[3];
rz(-1.3941744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.589754) q[2];
sx q[2];
rz(-1.7244312) q[2];
sx q[2];
rz(-3.0736249) q[2];
rz(3.0170278) q[3];
sx q[3];
rz(-0.3228651) q[3];
sx q[3];
rz(-1.386806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2215866) q[0];
sx q[0];
rz(-3.0060372) q[0];
sx q[0];
rz(2.8979229) q[0];
rz(-2.5098353) q[1];
sx q[1];
rz(-1.7383722) q[1];
sx q[1];
rz(1.3557281) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5403554) q[0];
sx q[0];
rz(-1.3155126) q[0];
sx q[0];
rz(3.113494) q[0];
rz(-0.83265702) q[2];
sx q[2];
rz(-0.50561935) q[2];
sx q[2];
rz(0.51072272) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6527378) q[1];
sx q[1];
rz(-1.6994209) q[1];
sx q[1];
rz(1.1210404) q[1];
x q[2];
rz(-2.5494266) q[3];
sx q[3];
rz(-1.6093996) q[3];
sx q[3];
rz(-1.5241227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0791066) q[2];
sx q[2];
rz(-2.1495543) q[2];
sx q[2];
rz(-0.24965723) q[2];
rz(2.6349973) q[3];
sx q[3];
rz(-1.6258312) q[3];
sx q[3];
rz(-0.33199582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8963985) q[0];
sx q[0];
rz(-1.2250552) q[0];
sx q[0];
rz(0.89843345) q[0];
rz(1.3348745) q[1];
sx q[1];
rz(-1.9060262) q[1];
sx q[1];
rz(1.2737087) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35667426) q[0];
sx q[0];
rz(-1.9264364) q[0];
sx q[0];
rz(-0.26892923) q[0];
rz(-1.6689698) q[2];
sx q[2];
rz(-1.86491) q[2];
sx q[2];
rz(-1.3054747) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.85325235) q[1];
sx q[1];
rz(-0.77008343) q[1];
sx q[1];
rz(-2.6283162) q[1];
rz(-2.0687194) q[3];
sx q[3];
rz(-0.91091279) q[3];
sx q[3];
rz(-2.4592196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0818103) q[2];
sx q[2];
rz(-0.60478294) q[2];
sx q[2];
rz(-0.95345062) q[2];
rz(-0.034514286) q[3];
sx q[3];
rz(-0.78648609) q[3];
sx q[3];
rz(-0.22687337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8614486) q[0];
sx q[0];
rz(-0.21629688) q[0];
sx q[0];
rz(0.24818534) q[0];
rz(1.0379627) q[1];
sx q[1];
rz(-2.018441) q[1];
sx q[1];
rz(-0.074137069) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5384597) q[0];
sx q[0];
rz(-2.5939301) q[0];
sx q[0];
rz(2.0752226) q[0];
rz(-pi) q[1];
rz(2.2494227) q[2];
sx q[2];
rz(-1.2954419) q[2];
sx q[2];
rz(-2.8868669) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.46596913) q[1];
sx q[1];
rz(-0.33756653) q[1];
sx q[1];
rz(1.9338495) q[1];
x q[2];
rz(-0.62319237) q[3];
sx q[3];
rz(-0.29286256) q[3];
sx q[3];
rz(2.5588536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2512102) q[2];
sx q[2];
rz(-2.7373098) q[2];
sx q[2];
rz(-3.1029491) q[2];
rz(0.97366992) q[3];
sx q[3];
rz(-0.49574167) q[3];
sx q[3];
rz(-2.8715449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6073109) q[0];
sx q[0];
rz(-1.5357635) q[0];
sx q[0];
rz(1.3624396) q[0];
rz(-0.81659395) q[1];
sx q[1];
rz(-1.2885619) q[1];
sx q[1];
rz(-1.1626676) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8602596) q[0];
sx q[0];
rz(-1.5409924) q[0];
sx q[0];
rz(1.6822862) q[0];
rz(-1.0983724) q[2];
sx q[2];
rz(-0.35944164) q[2];
sx q[2];
rz(-3.0430832) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3711277) q[1];
sx q[1];
rz(-2.7264997) q[1];
sx q[1];
rz(1.0789372) q[1];
rz(-pi) q[2];
rz(-0.88697042) q[3];
sx q[3];
rz(-2.5081222) q[3];
sx q[3];
rz(0.57852832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6267307) q[2];
sx q[2];
rz(-2.0613487) q[2];
sx q[2];
rz(-2.999372) q[2];
rz(2.2375315) q[3];
sx q[3];
rz(-1.3198493) q[3];
sx q[3];
rz(0.18946762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11480039) q[0];
sx q[0];
rz(-0.27652201) q[0];
sx q[0];
rz(-1.6739155) q[0];
rz(-0.57178512) q[1];
sx q[1];
rz(-0.3586868) q[1];
sx q[1];
rz(-0.30803672) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15333262) q[0];
sx q[0];
rz(-1.6969661) q[0];
sx q[0];
rz(-1.6660965) q[0];
rz(-pi) q[1];
rz(0.61200895) q[2];
sx q[2];
rz(-2.0888121) q[2];
sx q[2];
rz(0.81105622) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6412515) q[1];
sx q[1];
rz(-2.2441823) q[1];
sx q[1];
rz(2.889269) q[1];
x q[2];
rz(-0.34835784) q[3];
sx q[3];
rz(-1.6568686) q[3];
sx q[3];
rz(2.0406046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3623111) q[2];
sx q[2];
rz(-1.7411391) q[2];
sx q[2];
rz(1.9936838) q[2];
rz(2.4273196) q[3];
sx q[3];
rz(-2.2439984) q[3];
sx q[3];
rz(-2.4961297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66184735) q[0];
sx q[0];
rz(-2.3006738) q[0];
sx q[0];
rz(3.0116144) q[0];
rz(-0.030844363) q[1];
sx q[1];
rz(-1.8519311) q[1];
sx q[1];
rz(2.470509) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6993461) q[0];
sx q[0];
rz(-1.5171577) q[0];
sx q[0];
rz(1.5357114) q[0];
rz(-pi) q[1];
rz(1.408377) q[2];
sx q[2];
rz(-0.186609) q[2];
sx q[2];
rz(-2.9090372) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.259686) q[1];
sx q[1];
rz(-2.7555008) q[1];
sx q[1];
rz(1.8465471) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2460849) q[3];
sx q[3];
rz(-1.607778) q[3];
sx q[3];
rz(2.0946338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.122763) q[2];
sx q[2];
rz(-1.5315703) q[2];
sx q[2];
rz(0.57787952) q[2];
rz(0.028586483) q[3];
sx q[3];
rz(-1.280602) q[3];
sx q[3];
rz(-1.8813429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9776483) q[0];
sx q[0];
rz(-0.90286911) q[0];
sx q[0];
rz(0.41241616) q[0];
rz(-1.4498129) q[1];
sx q[1];
rz(-1.7990566) q[1];
sx q[1];
rz(-1.1669881) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3827688) q[0];
sx q[0];
rz(-2.1817657) q[0];
sx q[0];
rz(0.2610892) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86161676) q[2];
sx q[2];
rz(-1.78252) q[2];
sx q[2];
rz(0.35702969) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1046909) q[1];
sx q[1];
rz(-0.19100405) q[1];
sx q[1];
rz(0.16049338) q[1];
rz(-0.1498296) q[3];
sx q[3];
rz(-1.0458046) q[3];
sx q[3];
rz(-0.52469745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2010487) q[2];
sx q[2];
rz(-0.87934914) q[2];
sx q[2];
rz(-2.9525625) q[2];
rz(0.14686251) q[3];
sx q[3];
rz(-2.9569914) q[3];
sx q[3];
rz(1.7485025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015633164) q[0];
sx q[0];
rz(-1.8305612) q[0];
sx q[0];
rz(-0.92700672) q[0];
rz(1.3828297) q[1];
sx q[1];
rz(-2.5320876) q[1];
sx q[1];
rz(1.4896726) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2531567) q[0];
sx q[0];
rz(-2.2373767) q[0];
sx q[0];
rz(0.96320926) q[0];
rz(-1.6696879) q[2];
sx q[2];
rz(-2.6675468) q[2];
sx q[2];
rz(0.74741077) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5786963) q[1];
sx q[1];
rz(-1.3778731) q[1];
sx q[1];
rz(-1.3955411) q[1];
rz(-pi) q[2];
rz(-0.82724039) q[3];
sx q[3];
rz(-2.4628371) q[3];
sx q[3];
rz(0.37952207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5513409) q[2];
sx q[2];
rz(-0.44289032) q[2];
sx q[2];
rz(-1.4302953) q[2];
rz(0.57724214) q[3];
sx q[3];
rz(-0.8876628) q[3];
sx q[3];
rz(1.757471) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5532613) q[0];
sx q[0];
rz(-1.7734779) q[0];
sx q[0];
rz(-2.8531895) q[0];
rz(0.53238955) q[1];
sx q[1];
rz(-0.45982292) q[1];
sx q[1];
rz(-2.9945701) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64919103) q[0];
sx q[0];
rz(-2.1426139) q[0];
sx q[0];
rz(-0.4581106) q[0];
x q[1];
rz(-2.8617919) q[2];
sx q[2];
rz(-2.8923312) q[2];
sx q[2];
rz(-1.3485497) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.22132561) q[1];
sx q[1];
rz(-0.95278554) q[1];
sx q[1];
rz(-3.0926535) q[1];
rz(-pi) q[2];
rz(2.5649928) q[3];
sx q[3];
rz(-1.4398265) q[3];
sx q[3];
rz(3.1104345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0599351) q[2];
sx q[2];
rz(-0.43490484) q[2];
sx q[2];
rz(2.3975513) q[2];
rz(-2.384281) q[3];
sx q[3];
rz(-1.363874) q[3];
sx q[3];
rz(-1.7448759) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.025678) q[0];
sx q[0];
rz(-2.0712576) q[0];
sx q[0];
rz(2.0448137) q[0];
rz(2.3241282) q[1];
sx q[1];
rz(-1.2066963) q[1];
sx q[1];
rz(-0.6304601) q[1];
rz(3.0030737) q[2];
sx q[2];
rz(-0.45590966) q[2];
sx q[2];
rz(1.6675303) q[2];
rz(-1.7759454) q[3];
sx q[3];
rz(-2.57006) q[3];
sx q[3];
rz(-0.80001696) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
