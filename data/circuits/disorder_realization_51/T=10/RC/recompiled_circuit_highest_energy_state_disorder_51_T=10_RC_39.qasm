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
rz(-2.812204) q[0];
sx q[0];
rz(-0.63151276) q[0];
sx q[0];
rz(3.1407177) q[0];
rz(2.5007091) q[1];
sx q[1];
rz(-2.1407318) q[1];
sx q[1];
rz(-0.34520087) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5986818) q[0];
sx q[0];
rz(-0.094176725) q[0];
sx q[0];
rz(1.3046632) q[0];
x q[1];
rz(-0.35635524) q[2];
sx q[2];
rz(-0.40268597) q[2];
sx q[2];
rz(2.8860983) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.4243489) q[1];
sx q[1];
rz(-0.97704889) q[1];
sx q[1];
rz(0.51148606) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1253417) q[3];
sx q[3];
rz(-1.4964087) q[3];
sx q[3];
rz(-2.1824238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22902809) q[2];
sx q[2];
rz(-1.8077069) q[2];
sx q[2];
rz(2.933617) q[2];
rz(-0.29911706) q[3];
sx q[3];
rz(-2.5695473) q[3];
sx q[3];
rz(2.0007029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3308554) q[0];
sx q[0];
rz(-0.24092291) q[0];
sx q[0];
rz(-2.148707) q[0];
rz(-1.8244686) q[1];
sx q[1];
rz(-0.31610745) q[1];
sx q[1];
rz(2.3670926) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74906384) q[0];
sx q[0];
rz(-1.3925902) q[0];
sx q[0];
rz(0.011179608) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0335017) q[2];
sx q[2];
rz(-2.3024493) q[2];
sx q[2];
rz(3.0947859) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.61770536) q[1];
sx q[1];
rz(-1.7951709) q[1];
sx q[1];
rz(-0.38274204) q[1];
x q[2];
rz(-2.9162972) q[3];
sx q[3];
rz(-2.0301798) q[3];
sx q[3];
rz(-2.6090906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2449067) q[2];
sx q[2];
rz(-1.8550355) q[2];
sx q[2];
rz(-2.1265105) q[2];
rz(2.8924938) q[3];
sx q[3];
rz(-0.86779147) q[3];
sx q[3];
rz(2.7740313) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86431137) q[0];
sx q[0];
rz(-2.4503777) q[0];
sx q[0];
rz(-2.9840898) q[0];
rz(-1.0109673) q[1];
sx q[1];
rz(-0.33282655) q[1];
sx q[1];
rz(-3.0027622) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75051266) q[0];
sx q[0];
rz(-1.2475999) q[0];
sx q[0];
rz(-1.1884407) q[0];
rz(-pi) q[1];
rz(2.4262397) q[2];
sx q[2];
rz(-1.2497447) q[2];
sx q[2];
rz(-2.3356444) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0142747) q[1];
sx q[1];
rz(-1.7671314) q[1];
sx q[1];
rz(-1.2731958) q[1];
rz(-pi) q[2];
rz(1.3957455) q[3];
sx q[3];
rz(-1.850046) q[3];
sx q[3];
rz(1.4917013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.46131721) q[2];
sx q[2];
rz(-0.91247827) q[2];
sx q[2];
rz(0.3581363) q[2];
rz(2.4628468) q[3];
sx q[3];
rz(-2.1923784) q[3];
sx q[3];
rz(-2.1024735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.045227483) q[0];
sx q[0];
rz(-2.1684833) q[0];
sx q[0];
rz(-2.8367693) q[0];
rz(2.7335956) q[1];
sx q[1];
rz(-1.4381189) q[1];
sx q[1];
rz(-0.92794424) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.613425) q[0];
sx q[0];
rz(-1.6990464) q[0];
sx q[0];
rz(1.7665461) q[0];
rz(-2.1963652) q[2];
sx q[2];
rz(-1.06377) q[2];
sx q[2];
rz(2.8692109) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.18670652) q[1];
sx q[1];
rz(-0.88973239) q[1];
sx q[1];
rz(-1.7658556) q[1];
x q[2];
rz(-2.5119588) q[3];
sx q[3];
rz(-1.140652) q[3];
sx q[3];
rz(0.54508524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1912332) q[2];
sx q[2];
rz(-2.5706036) q[2];
sx q[2];
rz(0.18816571) q[2];
rz(-1.2763216) q[3];
sx q[3];
rz(-1.3151582) q[3];
sx q[3];
rz(1.7108542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-1.6998049) q[0];
sx q[0];
rz(-2.7975174) q[0];
sx q[0];
rz(-0.019388327) q[0];
rz(-1.060932) q[1];
sx q[1];
rz(-1.7273936) q[1];
sx q[1];
rz(-1.9020938) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0584377) q[0];
sx q[0];
rz(-2.4920336) q[0];
sx q[0];
rz(1.0970647) q[0];
rz(0.73915787) q[2];
sx q[2];
rz(-2.0817167) q[2];
sx q[2];
rz(-2.7834653) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6790948) q[1];
sx q[1];
rz(-1.0626918) q[1];
sx q[1];
rz(1.2485882) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89770384) q[3];
sx q[3];
rz(-2.5335059) q[3];
sx q[3];
rz(-1.9185818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1218607) q[2];
sx q[2];
rz(-1.3132361) q[2];
sx q[2];
rz(1.8149553) q[2];
rz(-0.2462247) q[3];
sx q[3];
rz(-2.0387869) q[3];
sx q[3];
rz(0.8518014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6061848) q[0];
sx q[0];
rz(-2.1562205) q[0];
sx q[0];
rz(-0.73915172) q[0];
rz(-0.34573653) q[1];
sx q[1];
rz(-1.3290936) q[1];
sx q[1];
rz(0.87127042) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3987253) q[0];
sx q[0];
rz(-2.350432) q[0];
sx q[0];
rz(-1.6683116) q[0];
rz(-2.3910693) q[2];
sx q[2];
rz(-2.0156246) q[2];
sx q[2];
rz(-1.9409279) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4193486) q[1];
sx q[1];
rz(-2.7461395) q[1];
sx q[1];
rz(-1.3802746) q[1];
rz(-pi) q[2];
rz(-2.0454117) q[3];
sx q[3];
rz(-2.1283009) q[3];
sx q[3];
rz(-1.6826009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8295916) q[2];
sx q[2];
rz(-0.63903725) q[2];
sx q[2];
rz(-2.269022) q[2];
rz(-1.568659) q[3];
sx q[3];
rz(-0.60397732) q[3];
sx q[3];
rz(0.93305552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.0503814) q[0];
sx q[0];
rz(-2.8822883) q[0];
sx q[0];
rz(-2.3799489) q[0];
rz(2.7161982) q[1];
sx q[1];
rz(-1.1401221) q[1];
sx q[1];
rz(2.2147307) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.923553) q[0];
sx q[0];
rz(-1.3430376) q[0];
sx q[0];
rz(-2.8995598) q[0];
x q[1];
rz(-0.71452801) q[2];
sx q[2];
rz(-2.0863669) q[2];
sx q[2];
rz(0.58802468) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4372448) q[1];
sx q[1];
rz(-1.3489375) q[1];
sx q[1];
rz(-1.941844) q[1];
x q[2];
rz(-2.1939932) q[3];
sx q[3];
rz(-1.5316846) q[3];
sx q[3];
rz(0.87621237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0285792) q[2];
sx q[2];
rz(-0.69602746) q[2];
sx q[2];
rz(-0.87116233) q[2];
rz(0.29843676) q[3];
sx q[3];
rz(-1.2454183) q[3];
sx q[3];
rz(1.3309853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7262064) q[0];
sx q[0];
rz(-2.6415249) q[0];
sx q[0];
rz(0.4183847) q[0];
rz(-2.8437974) q[1];
sx q[1];
rz(-1.1957542) q[1];
sx q[1];
rz(0.13430886) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1537841) q[0];
sx q[0];
rz(-1.0414818) q[0];
sx q[0];
rz(2.3321926) q[0];
rz(0.18608002) q[2];
sx q[2];
rz(-1.1717516) q[2];
sx q[2];
rz(0.068313561) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5843552) q[1];
sx q[1];
rz(-1.3882625) q[1];
sx q[1];
rz(1.7076769) q[1];
rz(-1.1014492) q[3];
sx q[3];
rz(-1.3162287) q[3];
sx q[3];
rz(-2.1397487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.40992752) q[2];
sx q[2];
rz(-1.2890559) q[2];
sx q[2];
rz(-2.4995787) q[2];
rz(-1.0632473) q[3];
sx q[3];
rz(-0.65170538) q[3];
sx q[3];
rz(-1.0412019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6006271) q[0];
sx q[0];
rz(-3.0525115) q[0];
sx q[0];
rz(-3.0294321) q[0];
rz(-1.8070096) q[1];
sx q[1];
rz(-0.59252512) q[1];
sx q[1];
rz(2.1122011) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5350069) q[0];
sx q[0];
rz(-0.86678737) q[0];
sx q[0];
rz(3.0993942) q[0];
rz(3.0156187) q[2];
sx q[2];
rz(-2.8354037) q[2];
sx q[2];
rz(0.15737113) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9255213) q[1];
sx q[1];
rz(-0.23434429) q[1];
sx q[1];
rz(1.3804803) q[1];
rz(-pi) q[2];
rz(0.46799) q[3];
sx q[3];
rz(-1.6525998) q[3];
sx q[3];
rz(1.8379267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42334291) q[2];
sx q[2];
rz(-0.074904718) q[2];
sx q[2];
rz(3.1035799) q[2];
rz(0.612261) q[3];
sx q[3];
rz(-2.2050048) q[3];
sx q[3];
rz(3.0003701) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8925979) q[0];
sx q[0];
rz(-1.4709512) q[0];
sx q[0];
rz(-0.41879642) q[0];
rz(-1.2576125) q[1];
sx q[1];
rz(-1.5092756) q[1];
sx q[1];
rz(-0.14258252) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2217956) q[0];
sx q[0];
rz(-0.16079535) q[0];
sx q[0];
rz(-0.17531403) q[0];
rz(1.293574) q[2];
sx q[2];
rz(-0.96053329) q[2];
sx q[2];
rz(2.0137613) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.31227641) q[1];
sx q[1];
rz(-1.8354776) q[1];
sx q[1];
rz(2.1618202) q[1];
rz(-pi) q[2];
rz(2.1964588) q[3];
sx q[3];
rz(-2.5925853) q[3];
sx q[3];
rz(1.4956724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7964145) q[2];
sx q[2];
rz(-3.0609481) q[2];
sx q[2];
rz(-2.8734015) q[2];
rz(-1.4043407) q[3];
sx q[3];
rz(-2.0495448) q[3];
sx q[3];
rz(-2.0973189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0161229) q[0];
sx q[0];
rz(-2.7653427) q[0];
sx q[0];
rz(0.34633037) q[0];
rz(-3.012433) q[1];
sx q[1];
rz(-1.2551413) q[1];
sx q[1];
rz(-1.7358949) q[1];
rz(-0.62025537) q[2];
sx q[2];
rz(-2.5778985) q[2];
sx q[2];
rz(2.3621205) q[2];
rz(3.075243) q[3];
sx q[3];
rz(-1.8649615) q[3];
sx q[3];
rz(-2.1129114) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
