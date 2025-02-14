OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2607245) q[0];
sx q[0];
rz(-0.27624929) q[0];
sx q[0];
rz(-2.2978388) q[0];
rz(-2.1387956) q[1];
sx q[1];
rz(-0.83477867) q[1];
sx q[1];
rz(-0.53258449) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0898429) q[0];
sx q[0];
rz(-2.2179363) q[0];
sx q[0];
rz(-0.18741742) q[0];
x q[1];
rz(-1.7108828) q[2];
sx q[2];
rz(-1.0647961) q[2];
sx q[2];
rz(-1.2337453) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4733359) q[1];
sx q[1];
rz(-1.419626) q[1];
sx q[1];
rz(0.87398793) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79382503) q[3];
sx q[3];
rz(-0.50758368) q[3];
sx q[3];
rz(2.6519597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8325995) q[2];
sx q[2];
rz(-1.4048046) q[2];
sx q[2];
rz(-0.12283202) q[2];
rz(3.099856) q[3];
sx q[3];
rz(-1.5580956) q[3];
sx q[3];
rz(-2.6532555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-1.1500583) q[0];
sx q[0];
rz(-2.8128862) q[0];
sx q[0];
rz(-1.1154037) q[0];
rz(-0.54488048) q[1];
sx q[1];
rz(-1.3976401) q[1];
sx q[1];
rz(0.56745183) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3504336) q[0];
sx q[0];
rz(-1.5251524) q[0];
sx q[0];
rz(1.6045536) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1092875) q[2];
sx q[2];
rz(-2.3650993) q[2];
sx q[2];
rz(2.1510498) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8711155) q[1];
sx q[1];
rz(-1.0024299) q[1];
sx q[1];
rz(1.4517484) q[1];
rz(0.53506923) q[3];
sx q[3];
rz(-2.5949083) q[3];
sx q[3];
rz(-3.1233623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6313717) q[2];
sx q[2];
rz(-3.0192182) q[2];
sx q[2];
rz(-3.0369634) q[2];
rz(2.6369324) q[3];
sx q[3];
rz(-1.3480836) q[3];
sx q[3];
rz(-2.4095355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.0013244) q[0];
sx q[0];
rz(-0.19135419) q[0];
sx q[0];
rz(-1.485317) q[0];
rz(-2.5775919) q[1];
sx q[1];
rz(-1.0537078) q[1];
sx q[1];
rz(2.3174813) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1230319) q[0];
sx q[0];
rz(-1.0433974) q[0];
sx q[0];
rz(-0.19077905) q[0];
rz(-pi) q[1];
rz(0.79470174) q[2];
sx q[2];
rz(-1.1247171) q[2];
sx q[2];
rz(2.5925555) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4802758) q[1];
sx q[1];
rz(-1.3696028) q[1];
sx q[1];
rz(-2.7072705) q[1];
rz(1.885627) q[3];
sx q[3];
rz(-2.9187429) q[3];
sx q[3];
rz(2.1999782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.04909) q[2];
sx q[2];
rz(-2.8503032) q[2];
sx q[2];
rz(1.606344) q[2];
rz(0.90318471) q[3];
sx q[3];
rz(-1.817037) q[3];
sx q[3];
rz(-2.7170392) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20151888) q[0];
sx q[0];
rz(-2.1551082) q[0];
sx q[0];
rz(0.41341138) q[0];
rz(-2.9460733) q[1];
sx q[1];
rz(-2.1045411) q[1];
sx q[1];
rz(1.6874541) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2233596) q[0];
sx q[0];
rz(-0.70230316) q[0];
sx q[0];
rz(-1.3732409) q[0];
rz(-pi) q[1];
rz(-0.35572534) q[2];
sx q[2];
rz(-1.7497831) q[2];
sx q[2];
rz(-2.9236874) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.153943) q[1];
sx q[1];
rz(-2.5763303) q[1];
sx q[1];
rz(-1.0260737) q[1];
x q[2];
rz(0.24202427) q[3];
sx q[3];
rz(-2.4325662) q[3];
sx q[3];
rz(0.55870648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.89895189) q[2];
sx q[2];
rz(-0.61598888) q[2];
sx q[2];
rz(-1.2658489) q[2];
rz(0.85593456) q[3];
sx q[3];
rz(-1.3866235) q[3];
sx q[3];
rz(1.9267193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6619381) q[0];
sx q[0];
rz(-2.9582773) q[0];
sx q[0];
rz(1.4187752) q[0];
rz(1.7105557) q[1];
sx q[1];
rz(-1.6173671) q[1];
sx q[1];
rz(0.72351825) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54004266) q[0];
sx q[0];
rz(-1.3941843) q[0];
sx q[0];
rz(0.58797449) q[0];
rz(-2.3846912) q[2];
sx q[2];
rz(-2.160103) q[2];
sx q[2];
rz(0.85739955) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0192467) q[1];
sx q[1];
rz(-0.86402383) q[1];
sx q[1];
rz(3.138423) q[1];
rz(-pi) q[2];
rz(1.5780903) q[3];
sx q[3];
rz(-1.2185343) q[3];
sx q[3];
rz(-1.1895113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7001223) q[2];
sx q[2];
rz(-1.0157601) q[2];
sx q[2];
rz(-0.36386841) q[2];
rz(-1.3077334) q[3];
sx q[3];
rz(-1.6677083) q[3];
sx q[3];
rz(-1.3522805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7135007) q[0];
sx q[0];
rz(-2.2245753) q[0];
sx q[0];
rz(-0.8999024) q[0];
rz(-2.0948386) q[1];
sx q[1];
rz(-2.7622107) q[1];
sx q[1];
rz(2.5894763) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13487694) q[0];
sx q[0];
rz(-1.9089886) q[0];
sx q[0];
rz(1.6007559) q[0];
x q[1];
rz(1.029676) q[2];
sx q[2];
rz(-2.1976802) q[2];
sx q[2];
rz(2.1205001) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.49657288) q[1];
sx q[1];
rz(-1.3995808) q[1];
sx q[1];
rz(0.091793072) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5512929) q[3];
sx q[3];
rz(-1.5029782) q[3];
sx q[3];
rz(2.6026311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.8206574) q[2];
sx q[2];
rz(-0.67639095) q[2];
sx q[2];
rz(-2.9429842) q[2];
rz(1.1292388) q[3];
sx q[3];
rz(-1.6695453) q[3];
sx q[3];
rz(-2.3614597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(2.8010913) q[0];
sx q[0];
rz(-1.3022364) q[0];
sx q[0];
rz(2.4275725) q[0];
rz(1.2227614) q[1];
sx q[1];
rz(-2.265265) q[1];
sx q[1];
rz(-2.8661935) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31322436) q[0];
sx q[0];
rz(-1.3592915) q[0];
sx q[0];
rz(1.9755892) q[0];
x q[1];
rz(-2.0049663) q[2];
sx q[2];
rz(-2.482802) q[2];
sx q[2];
rz(-0.24497797) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.95412817) q[1];
sx q[1];
rz(-0.70425941) q[1];
sx q[1];
rz(0.7974979) q[1];
x q[2];
rz(2.4163298) q[3];
sx q[3];
rz(-1.5607087) q[3];
sx q[3];
rz(-0.46848224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8347281) q[2];
sx q[2];
rz(-1.4707668) q[2];
sx q[2];
rz(-1.6738221) q[2];
rz(1.62014) q[3];
sx q[3];
rz(-0.68632564) q[3];
sx q[3];
rz(-2.2036688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15243212) q[0];
sx q[0];
rz(-2.2280362) q[0];
sx q[0];
rz(2.1754225) q[0];
rz(-2.3518111) q[1];
sx q[1];
rz(-2.8826394) q[1];
sx q[1];
rz(-0.97377473) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6286929) q[0];
sx q[0];
rz(-1.8159465) q[0];
sx q[0];
rz(-0.12314491) q[0];
x q[1];
rz(0.42736407) q[2];
sx q[2];
rz(-1.2334497) q[2];
sx q[2];
rz(0.34572476) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5681947) q[1];
sx q[1];
rz(-0.97028448) q[1];
sx q[1];
rz(-2.7706514) q[1];
rz(1.4548746) q[3];
sx q[3];
rz(-0.478607) q[3];
sx q[3];
rz(-2.6461864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8345653) q[2];
sx q[2];
rz(-1.4296738) q[2];
sx q[2];
rz(-0.60565051) q[2];
rz(-1.0904795) q[3];
sx q[3];
rz(-0.24221371) q[3];
sx q[3];
rz(-3.0428913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.403648) q[0];
sx q[0];
rz(-3.0927959) q[0];
sx q[0];
rz(0.57156372) q[0];
rz(-0.38772186) q[1];
sx q[1];
rz(-0.78920904) q[1];
sx q[1];
rz(-0.8943843) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1739367) q[0];
sx q[0];
rz(-1.0595788) q[0];
sx q[0];
rz(2.5777667) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11668514) q[2];
sx q[2];
rz(-2.1348079) q[2];
sx q[2];
rz(2.2224768) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5928985) q[1];
sx q[1];
rz(-1.1709524) q[1];
sx q[1];
rz(-0.26229761) q[1];
x q[2];
rz(-0.13988755) q[3];
sx q[3];
rz(-2.6287946) q[3];
sx q[3];
rz(-0.96641216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9992708) q[2];
sx q[2];
rz(-0.88185507) q[2];
sx q[2];
rz(-2.2692661) q[2];
rz(-0.074660389) q[3];
sx q[3];
rz(-1.9118237) q[3];
sx q[3];
rz(-2.5653896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9999303) q[0];
sx q[0];
rz(-0.44121989) q[0];
sx q[0];
rz(-2.4040851) q[0];
rz(-2.4420338) q[1];
sx q[1];
rz(-1.0195426) q[1];
sx q[1];
rz(-0.84651822) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2071211) q[0];
sx q[0];
rz(-1.9686598) q[0];
sx q[0];
rz(-0.69402309) q[0];
rz(-pi) q[1];
rz(-3.0902683) q[2];
sx q[2];
rz(-1.996417) q[2];
sx q[2];
rz(-1.0251306) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6168931) q[1];
sx q[1];
rz(-1.9977101) q[1];
sx q[1];
rz(-2.9906143) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3844526) q[3];
sx q[3];
rz(-1.1769503) q[3];
sx q[3];
rz(-2.9316616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6503341) q[2];
sx q[2];
rz(-1.1639872) q[2];
sx q[2];
rz(0.35438806) q[2];
rz(-2.3645511) q[3];
sx q[3];
rz(-0.6207501) q[3];
sx q[3];
rz(-2.0508155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1914094) q[0];
sx q[0];
rz(-1.8671028) q[0];
sx q[0];
rz(-2.6381459) q[0];
rz(1.3632111) q[1];
sx q[1];
rz(-2.6111205) q[1];
sx q[1];
rz(3.0725239) q[1];
rz(1.4466046) q[2];
sx q[2];
rz(-0.71659452) q[2];
sx q[2];
rz(2.334395) q[2];
rz(1.1152399) q[3];
sx q[3];
rz(-1.8805877) q[3];
sx q[3];
rz(2.7177802) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
