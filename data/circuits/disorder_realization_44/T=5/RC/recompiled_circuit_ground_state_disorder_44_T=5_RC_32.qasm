OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3888336) q[0];
sx q[0];
rz(-1.692481) q[0];
sx q[0];
rz(-1.4504855) q[0];
rz(0.9736355) q[1];
sx q[1];
rz(-1.7042301) q[1];
sx q[1];
rz(-0.91926423) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5333119) q[0];
sx q[0];
rz(-1.5402147) q[0];
sx q[0];
rz(-1.5536867) q[0];
rz(-3.0212901) q[2];
sx q[2];
rz(-1.4805111) q[2];
sx q[2];
rz(1.1995969) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3340993) q[1];
sx q[1];
rz(-1.6624644) q[1];
sx q[1];
rz(-2.2794072) q[1];
rz(-pi) q[2];
rz(1.1000865) q[3];
sx q[3];
rz(-1.9506427) q[3];
sx q[3];
rz(2.0947411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8093402) q[2];
sx q[2];
rz(-1.6925749) q[2];
sx q[2];
rz(-1.7285041) q[2];
rz(-0.20279065) q[3];
sx q[3];
rz(-1.3821802) q[3];
sx q[3];
rz(0.10281674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8973812) q[0];
sx q[0];
rz(-2.057071) q[0];
sx q[0];
rz(-2.4308423) q[0];
rz(-0.76849014) q[1];
sx q[1];
rz(-2.0748506) q[1];
sx q[1];
rz(-2.1220727) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43895129) q[0];
sx q[0];
rz(-1.1239237) q[0];
sx q[0];
rz(1.5269482) q[0];
x q[1];
rz(1.1395733) q[2];
sx q[2];
rz(-1.5016455) q[2];
sx q[2];
rz(0.56996843) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.72538917) q[1];
sx q[1];
rz(-0.37230154) q[1];
sx q[1];
rz(-1.4854234) q[1];
rz(-1.4367661) q[3];
sx q[3];
rz(-0.55826) q[3];
sx q[3];
rz(0.2414862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3780313) q[2];
sx q[2];
rz(-1.8969994) q[2];
sx q[2];
rz(1.2223318) q[2];
rz(-1.2126806) q[3];
sx q[3];
rz(-2.1247532) q[3];
sx q[3];
rz(2.4299664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0843622) q[0];
sx q[0];
rz(-0.98452345) q[0];
sx q[0];
rz(-1.2179751) q[0];
rz(0.51586622) q[1];
sx q[1];
rz(-0.54793826) q[1];
sx q[1];
rz(-0.9224433) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9003214) q[0];
sx q[0];
rz(-1.5301955) q[0];
sx q[0];
rz(-3.0837644) q[0];
rz(2.416196) q[2];
sx q[2];
rz(-1.179152) q[2];
sx q[2];
rz(1.894941) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6092012) q[1];
sx q[1];
rz(-0.89840404) q[1];
sx q[1];
rz(2.3228541) q[1];
x q[2];
rz(2.8124468) q[3];
sx q[3];
rz(-1.7308047) q[3];
sx q[3];
rz(-3.126006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1232274) q[2];
sx q[2];
rz(-2.0570698) q[2];
sx q[2];
rz(-1.8017192) q[2];
rz(-0.54723048) q[3];
sx q[3];
rz(-2.7180143) q[3];
sx q[3];
rz(-0.71119285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.055534) q[0];
sx q[0];
rz(-0.14500293) q[0];
sx q[0];
rz(-2.9823629) q[0];
rz(3.1314462) q[1];
sx q[1];
rz(-2.1187014) q[1];
sx q[1];
rz(-0.25746447) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5978569) q[0];
sx q[0];
rz(-2.275995) q[0];
sx q[0];
rz(-2.0446834) q[0];
rz(-pi) q[1];
x q[1];
rz(0.065071062) q[2];
sx q[2];
rz(-2.3461968) q[2];
sx q[2];
rz(-2.3415023) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3066669) q[1];
sx q[1];
rz(-1.67027) q[1];
sx q[1];
rz(2.5132781) q[1];
rz(0.61473989) q[3];
sx q[3];
rz(-0.90723824) q[3];
sx q[3];
rz(3.1082982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3782392) q[2];
sx q[2];
rz(-1.2946318) q[2];
sx q[2];
rz(-0.47284687) q[2];
rz(-0.7044479) q[3];
sx q[3];
rz(-1.7770146) q[3];
sx q[3];
rz(-0.64594597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6351629) q[0];
sx q[0];
rz(-1.9671054) q[0];
sx q[0];
rz(0.065486431) q[0];
rz(2.7323515) q[1];
sx q[1];
rz(-1.1301273) q[1];
sx q[1];
rz(1.4261036) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.586389) q[0];
sx q[0];
rz(-2.9034333) q[0];
sx q[0];
rz(-2.737153) q[0];
rz(-pi) q[1];
rz(-2.7528115) q[2];
sx q[2];
rz(-1.4884623) q[2];
sx q[2];
rz(1.3363234) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3196484) q[1];
sx q[1];
rz(-1.2927755) q[1];
sx q[1];
rz(-1.8026428) q[1];
rz(-pi) q[2];
rz(-2.0992005) q[3];
sx q[3];
rz(-2.8082153) q[3];
sx q[3];
rz(2.7957145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9466729) q[2];
sx q[2];
rz(-2.0544923) q[2];
sx q[2];
rz(-1.3067513) q[2];
rz(1.170018) q[3];
sx q[3];
rz(-2.7368059) q[3];
sx q[3];
rz(-0.10725966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4829247) q[0];
sx q[0];
rz(-1.1558477) q[0];
sx q[0];
rz(0.40147993) q[0];
rz(0.67277706) q[1];
sx q[1];
rz(-2.3428226) q[1];
sx q[1];
rz(2.2875517) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5082762) q[0];
sx q[0];
rz(-0.22249732) q[0];
sx q[0];
rz(1.7822687) q[0];
x q[1];
rz(-2.4638205) q[2];
sx q[2];
rz(-2.9052264) q[2];
sx q[2];
rz(1.7198228) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5062949) q[1];
sx q[1];
rz(-2.1599401) q[1];
sx q[1];
rz(-1.4511257) q[1];
x q[2];
rz(-2.5402732) q[3];
sx q[3];
rz(-0.70402217) q[3];
sx q[3];
rz(-2.6773334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4868769) q[2];
sx q[2];
rz(-2.2701023) q[2];
sx q[2];
rz(1.9640478) q[2];
rz(0.55142895) q[3];
sx q[3];
rz(-1.5827554) q[3];
sx q[3];
rz(-2.3059755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7198782) q[0];
sx q[0];
rz(-2.8604909) q[0];
sx q[0];
rz(1.9792492) q[0];
rz(-1.989919) q[1];
sx q[1];
rz(-1.3538227) q[1];
sx q[1];
rz(0.91845671) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29354039) q[0];
sx q[0];
rz(-1.9013202) q[0];
sx q[0];
rz(0.28384112) q[0];
x q[1];
rz(2.8807441) q[2];
sx q[2];
rz(-2.2493304) q[2];
sx q[2];
rz(-1.4229753) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9388401) q[1];
sx q[1];
rz(-1.5781286) q[1];
sx q[1];
rz(-1.5640902) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1226467) q[3];
sx q[3];
rz(-0.3780685) q[3];
sx q[3];
rz(2.700875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.034417001) q[2];
sx q[2];
rz(-1.7566661) q[2];
sx q[2];
rz(-2.6386063) q[2];
rz(-2.3668187) q[3];
sx q[3];
rz(-0.46190244) q[3];
sx q[3];
rz(1.3431965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4485432) q[0];
sx q[0];
rz(-0.22668426) q[0];
sx q[0];
rz(-1.5013129) q[0];
rz(0.34128183) q[1];
sx q[1];
rz(-1.4309859) q[1];
sx q[1];
rz(-2.0223845) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2687896) q[0];
sx q[0];
rz(-0.69376341) q[0];
sx q[0];
rz(2.540178) q[0];
rz(-1.0897175) q[2];
sx q[2];
rz(-1.8029034) q[2];
sx q[2];
rz(0.26929917) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.95204845) q[1];
sx q[1];
rz(-0.91390007) q[1];
sx q[1];
rz(-1.4178965) q[1];
x q[2];
rz(-2.8265727) q[3];
sx q[3];
rz(-2.7286988) q[3];
sx q[3];
rz(1.1668432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1559653) q[2];
sx q[2];
rz(-2.1817744) q[2];
sx q[2];
rz(0.36925527) q[2];
rz(2.8333832) q[3];
sx q[3];
rz(-0.9674558) q[3];
sx q[3];
rz(-1.5306028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8592598) q[0];
sx q[0];
rz(-1.3066602) q[0];
sx q[0];
rz(0.34307137) q[0];
rz(-1.9725017) q[1];
sx q[1];
rz(-1.7770551) q[1];
sx q[1];
rz(1.4287359) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0562623) q[0];
sx q[0];
rz(-1.82162) q[0];
sx q[0];
rz(-0.40928264) q[0];
x q[1];
rz(2.3071204) q[2];
sx q[2];
rz(-1.021046) q[2];
sx q[2];
rz(2.2817734) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8642042) q[1];
sx q[1];
rz(-0.38071796) q[1];
sx q[1];
rz(-1.7362795) q[1];
rz(0.87782209) q[3];
sx q[3];
rz(-1.9435427) q[3];
sx q[3];
rz(-2.8543548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2009361) q[2];
sx q[2];
rz(-2.9328465) q[2];
sx q[2];
rz(-1.7500056) q[2];
rz(0.40677795) q[3];
sx q[3];
rz(-1.6986366) q[3];
sx q[3];
rz(0.87882915) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10130356) q[0];
sx q[0];
rz(-0.60281301) q[0];
sx q[0];
rz(-1.783675) q[0];
rz(1.2376002) q[1];
sx q[1];
rz(-1.0126746) q[1];
sx q[1];
rz(-0.79992574) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3615204) q[0];
sx q[0];
rz(-1.8006386) q[0];
sx q[0];
rz(-0.80432463) q[0];
rz(-1.2679891) q[2];
sx q[2];
rz(-1.0410415) q[2];
sx q[2];
rz(-0.15632665) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.057803072) q[1];
sx q[1];
rz(-2.4531595) q[1];
sx q[1];
rz(2.8958984) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0217729) q[3];
sx q[3];
rz(-2.6116237) q[3];
sx q[3];
rz(-1.8161023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7633729) q[2];
sx q[2];
rz(-1.4241445) q[2];
sx q[2];
rz(3.0685032) q[2];
rz(0.25589219) q[3];
sx q[3];
rz(-0.81232324) q[3];
sx q[3];
rz(-1.6528486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8681317) q[0];
sx q[0];
rz(-1.68597) q[0];
sx q[0];
rz(-1.2706533) q[0];
rz(-1.801626) q[1];
sx q[1];
rz(-1.3506964) q[1];
sx q[1];
rz(0.66257308) q[1];
rz(-2.432178) q[2];
sx q[2];
rz(-1.5723036) q[2];
sx q[2];
rz(0.2607762) q[2];
rz(-2.5005199) q[3];
sx q[3];
rz(-2.0430293) q[3];
sx q[3];
rz(0.52386491) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
