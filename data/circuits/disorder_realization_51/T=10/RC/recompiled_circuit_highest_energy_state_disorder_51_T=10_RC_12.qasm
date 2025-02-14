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
rz(-0.00087498571) q[0];
rz(-0.64088351) q[1];
sx q[1];
rz(-1.0008608) q[1];
sx q[1];
rz(0.34520087) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54291081) q[0];
sx q[0];
rz(-3.0474159) q[0];
sx q[0];
rz(1.3046632) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35635524) q[2];
sx q[2];
rz(-0.40268597) q[2];
sx q[2];
rz(-0.25549437) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6908603) q[1];
sx q[1];
rz(-1.1530515) q[1];
sx q[1];
rz(-0.91207204) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4963989) q[3];
sx q[3];
rz(-1.5545903) q[3];
sx q[3];
rz(0.61283535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9125646) q[2];
sx q[2];
rz(-1.3338858) q[2];
sx q[2];
rz(-0.20797569) q[2];
rz(0.29911706) q[3];
sx q[3];
rz(-2.5695473) q[3];
sx q[3];
rz(1.1408898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8107373) q[0];
sx q[0];
rz(-2.9006697) q[0];
sx q[0];
rz(0.99288565) q[0];
rz(-1.8244686) q[1];
sx q[1];
rz(-2.8254852) q[1];
sx q[1];
rz(-2.3670926) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82371432) q[0];
sx q[0];
rz(-1.5817989) q[0];
sx q[0];
rz(1.3925793) q[0];
rz(-pi) q[1];
rz(-0.80747898) q[2];
sx q[2];
rz(-1.9614855) q[2];
sx q[2];
rz(1.9964068) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6933813) q[1];
sx q[1];
rz(-0.44084545) q[1];
sx q[1];
rz(-0.54852672) q[1];
x q[2];
rz(-1.1011613) q[3];
sx q[3];
rz(-1.3691994) q[3];
sx q[3];
rz(-0.93702173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2449067) q[2];
sx q[2];
rz(-1.2865571) q[2];
sx q[2];
rz(-2.1265105) q[2];
rz(-0.24909881) q[3];
sx q[3];
rz(-0.86779147) q[3];
sx q[3];
rz(-0.36756137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2772813) q[0];
sx q[0];
rz(-0.69121498) q[0];
sx q[0];
rz(2.9840898) q[0];
rz(2.1306254) q[1];
sx q[1];
rz(-2.8087661) q[1];
sx q[1];
rz(3.0027622) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.39108) q[0];
sx q[0];
rz(-1.8939928) q[0];
sx q[0];
rz(-1.9531519) q[0];
x q[1];
rz(-1.9857668) q[2];
sx q[2];
rz(-2.2425644) q[2];
sx q[2];
rz(-0.49723724) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2646059) q[1];
sx q[1];
rz(-2.786676) q[1];
sx q[1];
rz(-0.97477977) q[1];
rz(-0.54581996) q[3];
sx q[3];
rz(-2.8132317) q[3];
sx q[3];
rz(-0.92121802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6802754) q[2];
sx q[2];
rz(-2.2291144) q[2];
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
rz(0.045227483) q[0];
sx q[0];
rz(-2.1684833) q[0];
sx q[0];
rz(-0.3048234) q[0];
rz(2.7335956) q[1];
sx q[1];
rz(-1.7034737) q[1];
sx q[1];
rz(0.92794424) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6154895) q[0];
sx q[0];
rz(-0.23357059) q[0];
sx q[0];
rz(0.98532565) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5408542) q[2];
sx q[2];
rz(-2.1081807) q[2];
sx q[2];
rz(-1.5058277) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9548861) q[1];
sx q[1];
rz(-0.88973239) q[1];
sx q[1];
rz(-1.7658556) q[1];
rz(-pi) q[2];
rz(0.62963387) q[3];
sx q[3];
rz(-1.140652) q[3];
sx q[3];
rz(-2.5965074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1912332) q[2];
sx q[2];
rz(-2.5706036) q[2];
sx q[2];
rz(-0.18816571) q[2];
rz(1.865271) q[3];
sx q[3];
rz(-1.8264344) q[3];
sx q[3];
rz(-1.7108542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6998049) q[0];
sx q[0];
rz(-0.34407523) q[0];
sx q[0];
rz(0.019388327) q[0];
rz(-1.060932) q[1];
sx q[1];
rz(-1.414199) q[1];
sx q[1];
rz(-1.2394989) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083154924) q[0];
sx q[0];
rz(-0.64955901) q[0];
sx q[0];
rz(-1.0970647) q[0];
rz(-2.4024348) q[2];
sx q[2];
rz(-2.0817167) q[2];
sx q[2];
rz(0.3581274) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6790948) q[1];
sx q[1];
rz(-1.0626918) q[1];
sx q[1];
rz(1.8930045) q[1];
rz(-2.7321759) q[3];
sx q[3];
rz(-1.1077266) q[3];
sx q[3];
rz(2.6894929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0197319) q[2];
sx q[2];
rz(-1.3132361) q[2];
sx q[2];
rz(1.8149553) q[2];
rz(-2.895368) q[3];
sx q[3];
rz(-2.0387869) q[3];
sx q[3];
rz(-0.8518014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5354079) q[0];
sx q[0];
rz(-2.1562205) q[0];
sx q[0];
rz(-0.73915172) q[0];
rz(-0.34573653) q[1];
sx q[1];
rz(-1.3290936) q[1];
sx q[1];
rz(-2.2703222) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5369878) q[0];
sx q[0];
rz(-0.78442998) q[0];
sx q[0];
rz(3.0434199) q[0];
rz(-0.75052336) q[2];
sx q[2];
rz(-1.125968) q[2];
sx q[2];
rz(-1.9409279) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.516219) q[1];
sx q[1];
rz(-1.1828848) q[1];
sx q[1];
rz(0.078887786) q[1];
rz(-pi) q[2];
rz(-2.0454117) q[3];
sx q[3];
rz(-1.0132917) q[3];
sx q[3];
rz(1.6826009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.31200108) q[2];
sx q[2];
rz(-0.63903725) q[2];
sx q[2];
rz(-0.87257067) q[2];
rz(-1.5729337) q[3];
sx q[3];
rz(-0.60397732) q[3];
sx q[3];
rz(2.2085371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0912112) q[0];
sx q[0];
rz(-0.25930431) q[0];
sx q[0];
rz(2.3799489) q[0];
rz(-2.7161982) q[1];
sx q[1];
rz(-1.1401221) q[1];
sx q[1];
rz(0.92686191) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0480014) q[0];
sx q[0];
rz(-0.33080745) q[0];
sx q[0];
rz(-2.3729416) q[0];
rz(-pi) q[1];
rz(0.71304597) q[2];
sx q[2];
rz(-2.2879061) q[2];
sx q[2];
rz(-2.6756659) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.35173479) q[1];
sx q[1];
rz(-0.42966336) q[1];
sx q[1];
rz(2.1273145) q[1];
x q[2];
rz(-0.94759946) q[3];
sx q[3];
rz(-1.5316846) q[3];
sx q[3];
rz(2.2653803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0285792) q[2];
sx q[2];
rz(-0.69602746) q[2];
sx q[2];
rz(2.2704303) q[2];
rz(-0.29843676) q[3];
sx q[3];
rz(-1.8961743) q[3];
sx q[3];
rz(-1.8106073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7262064) q[0];
sx q[0];
rz(-0.50006777) q[0];
sx q[0];
rz(-0.4183847) q[0];
rz(0.2977953) q[1];
sx q[1];
rz(-1.1957542) q[1];
sx q[1];
rz(-3.0072838) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1537841) q[0];
sx q[0];
rz(-2.1001108) q[0];
sx q[0];
rz(2.3321926) q[0];
rz(-pi) q[1];
rz(2.9555126) q[2];
sx q[2];
rz(-1.1717516) q[2];
sx q[2];
rz(-0.068313561) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2063576) q[1];
sx q[1];
rz(-2.9138953) q[1];
sx q[1];
rz(-0.63657324) q[1];
x q[2];
rz(-1.0487324) q[3];
sx q[3];
rz(-2.6122141) q[3];
sx q[3];
rz(0.10806187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.40992752) q[2];
sx q[2];
rz(-1.2890559) q[2];
sx q[2];
rz(2.4995787) q[2];
rz(1.0632473) q[3];
sx q[3];
rz(-0.65170538) q[3];
sx q[3];
rz(1.0412019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5409656) q[0];
sx q[0];
rz(-0.0890812) q[0];
sx q[0];
rz(0.11216057) q[0];
rz(1.3345831) q[1];
sx q[1];
rz(-0.59252512) q[1];
sx q[1];
rz(2.1122011) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5350069) q[0];
sx q[0];
rz(-0.86678737) q[0];
sx q[0];
rz(-0.042198472) q[0];
x q[1];
rz(-1.5310982) q[2];
sx q[2];
rz(-1.8744812) q[2];
sx q[2];
rz(0.025321753) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9721016) q[1];
sx q[1];
rz(-1.5268561) q[1];
sx q[1];
rz(1.8010587) q[1];
rz(-pi) q[2];
rz(1.4791895) q[3];
sx q[3];
rz(-1.1044958) q[3];
sx q[3];
rz(-0.3084076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7182497) q[2];
sx q[2];
rz(-3.0666879) q[2];
sx q[2];
rz(-0.038012803) q[2];
rz(-0.612261) q[3];
sx q[3];
rz(-0.93658787) q[3];
sx q[3];
rz(3.0003701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24899471) q[0];
sx q[0];
rz(-1.4709512) q[0];
sx q[0];
rz(-0.41879642) q[0];
rz(1.8839802) q[1];
sx q[1];
rz(-1.6323171) q[1];
sx q[1];
rz(0.14258252) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.044238) q[0];
sx q[0];
rz(-1.4124866) q[0];
sx q[0];
rz(1.5425141) q[0];
x q[1];
rz(-2.5129287) q[2];
sx q[2];
rz(-1.3445879) q[2];
sx q[2];
rz(-2.5369801) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6302418) q[1];
sx q[1];
rz(-2.5005184) q[1];
sx q[1];
rz(2.0235233) q[1];
rz(-2.0311277) q[3];
sx q[3];
rz(-1.2602196) q[3];
sx q[3];
rz(0.47720695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7964145) q[2];
sx q[2];
rz(-0.08064457) q[2];
sx q[2];
rz(-0.26819116) q[2];
rz(1.7372519) q[3];
sx q[3];
rz(-2.0495448) q[3];
sx q[3];
rz(1.0442737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1254697) q[0];
sx q[0];
rz(-2.7653427) q[0];
sx q[0];
rz(0.34633037) q[0];
rz(-0.12915962) q[1];
sx q[1];
rz(-1.8864514) q[1];
sx q[1];
rz(1.4056978) q[1];
rz(-2.5213373) q[2];
sx q[2];
rz(-0.56369416) q[2];
sx q[2];
rz(-0.77947215) q[2];
rz(-1.3553452) q[3];
sx q[3];
rz(-2.8402495) q[3];
sx q[3];
rz(0.80339669) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
