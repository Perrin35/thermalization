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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54291081) q[0];
sx q[0];
rz(-3.0474159) q[0];
sx q[0];
rz(-1.8369294) q[0];
x q[1];
rz(0.35635524) q[2];
sx q[2];
rz(-0.40268597) q[2];
sx q[2];
rz(0.25549437) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6908603) q[1];
sx q[1];
rz(-1.1530515) q[1];
sx q[1];
rz(-2.2295206) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4963989) q[3];
sx q[3];
rz(-1.5545903) q[3];
sx q[3];
rz(-2.5287573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.22902809) q[2];
sx q[2];
rz(-1.3338858) q[2];
sx q[2];
rz(-2.933617) q[2];
rz(0.29911706) q[3];
sx q[3];
rz(-0.57204539) q[3];
sx q[3];
rz(-1.1408898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-1.3308554) q[0];
sx q[0];
rz(-0.24092291) q[0];
sx q[0];
rz(2.148707) q[0];
rz(-1.8244686) q[1];
sx q[1];
rz(-2.8254852) q[1];
sx q[1];
rz(-2.3670926) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3925288) q[0];
sx q[0];
rz(-1.3925902) q[0];
sx q[0];
rz(0.011179608) q[0];
rz(-0.51807816) q[2];
sx q[2];
rz(-0.87730125) q[2];
sx q[2];
rz(-0.77502807) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6933813) q[1];
sx q[1];
rz(-2.7007472) q[1];
sx q[1];
rz(2.5930659) q[1];
rz(-pi) q[2];
rz(-1.1466188) q[3];
sx q[3];
rz(-0.50809233) q[3];
sx q[3];
rz(1.0095694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2449067) q[2];
sx q[2];
rz(-1.2865571) q[2];
sx q[2];
rz(-2.1265105) q[2];
rz(2.8924938) q[3];
sx q[3];
rz(-2.2738012) q[3];
sx q[3];
rz(-2.7740313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2772813) q[0];
sx q[0];
rz(-0.69121498) q[0];
sx q[0];
rz(-2.9840898) q[0];
rz(1.0109673) q[1];
sx q[1];
rz(-0.33282655) q[1];
sx q[1];
rz(-0.13883042) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1942753) q[0];
sx q[0];
rz(-1.9324158) q[0];
sx q[0];
rz(2.7951434) q[0];
rz(-1.1558258) q[2];
sx q[2];
rz(-0.89902821) q[2];
sx q[2];
rz(2.6443554) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.638354) q[1];
sx q[1];
rz(-1.279083) q[1];
sx q[1];
rz(2.9364763) q[1];
rz(-0.54581996) q[3];
sx q[3];
rz(-0.32836093) q[3];
sx q[3];
rz(0.92121802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6802754) q[2];
sx q[2];
rz(-2.2291144) q[2];
sx q[2];
rz(-0.3581363) q[2];
rz(-0.67874587) q[3];
sx q[3];
rz(-0.94921422) q[3];
sx q[3];
rz(2.1024735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-0.045227483) q[0];
sx q[0];
rz(-2.1684833) q[0];
sx q[0];
rz(0.3048234) q[0];
rz(0.40799704) q[1];
sx q[1];
rz(-1.4381189) q[1];
sx q[1];
rz(0.92794424) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1243195) q[0];
sx q[0];
rz(-1.3766748) q[0];
sx q[0];
rz(-0.13071901) q[0];
x q[1];
rz(-0.81176968) q[2];
sx q[2];
rz(-0.78321811) q[2];
sx q[2];
rz(-0.70657544) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2603307) q[1];
sx q[1];
rz(-1.4196383) q[1];
sx q[1];
rz(-0.69044729) q[1];
rz(-pi) q[2];
rz(2.4797012) q[3];
sx q[3];
rz(-2.3958979) q[3];
sx q[3];
rz(2.6357366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1912332) q[2];
sx q[2];
rz(-2.5706036) q[2];
sx q[2];
rz(0.18816571) q[2];
rz(-1.2763216) q[3];
sx q[3];
rz(-1.8264344) q[3];
sx q[3];
rz(1.4307384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4417878) q[0];
sx q[0];
rz(-0.34407523) q[0];
sx q[0];
rz(3.1222043) q[0];
rz(2.0806606) q[1];
sx q[1];
rz(-1.414199) q[1];
sx q[1];
rz(1.9020938) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65514046) q[0];
sx q[0];
rz(-1.0024655) q[0];
sx q[0];
rz(-2.8080432) q[0];
x q[1];
rz(0.6940191) q[2];
sx q[2];
rz(-2.2713714) q[2];
sx q[2];
rz(-1.7050336) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2693138) q[1];
sx q[1];
rz(-1.8511103) q[1];
sx q[1];
rz(-0.53086908) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7321759) q[3];
sx q[3];
rz(-2.0338661) q[3];
sx q[3];
rz(2.6894929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1218607) q[2];
sx q[2];
rz(-1.3132361) q[2];
sx q[2];
rz(1.3266374) q[2];
rz(0.2462247) q[3];
sx q[3];
rz(-2.0387869) q[3];
sx q[3];
rz(2.2897913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5354079) q[0];
sx q[0];
rz(-0.98537213) q[0];
sx q[0];
rz(-2.4024409) q[0];
rz(0.34573653) q[1];
sx q[1];
rz(-1.3290936) q[1];
sx q[1];
rz(-0.87127042) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3987253) q[0];
sx q[0];
rz(-0.79116066) q[0];
sx q[0];
rz(-1.4732811) q[0];
x q[1];
rz(0.99314697) q[2];
sx q[2];
rz(-2.2340074) q[2];
sx q[2];
rz(-3.1301067) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.024684357) q[1];
sx q[1];
rz(-1.4977807) q[1];
sx q[1];
rz(-1.9597998) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5302251) q[3];
sx q[3];
rz(-1.9690367) q[3];
sx q[3];
rz(-2.7643725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.31200108) q[2];
sx q[2];
rz(-0.63903725) q[2];
sx q[2];
rz(-2.269022) q[2];
rz(-1.5729337) q[3];
sx q[3];
rz(-0.60397732) q[3];
sx q[3];
rz(2.2085371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(2.0912112) q[0];
sx q[0];
rz(-2.8822883) q[0];
sx q[0];
rz(2.3799489) q[0];
rz(-0.42539445) q[1];
sx q[1];
rz(-1.1401221) q[1];
sx q[1];
rz(-0.92686191) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8445209) q[0];
sx q[0];
rz(-1.8064587) q[0];
sx q[0];
rz(-1.8051487) q[0];
rz(-2.4285467) q[2];
sx q[2];
rz(-0.8536866) q[2];
sx q[2];
rz(-0.46592679) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4372448) q[1];
sx q[1];
rz(-1.7926551) q[1];
sx q[1];
rz(-1.941844) q[1];
x q[2];
rz(-1.5038483) q[3];
sx q[3];
rz(-2.5173325) q[3];
sx q[3];
rz(2.3926596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1130134) q[2];
sx q[2];
rz(-2.4455652) q[2];
sx q[2];
rz(0.87116233) q[2];
rz(0.29843676) q[3];
sx q[3];
rz(-1.2454183) q[3];
sx q[3];
rz(-1.8106073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7262064) q[0];
sx q[0];
rz(-0.50006777) q[0];
sx q[0];
rz(-0.4183847) q[0];
rz(-0.2977953) q[1];
sx q[1];
rz(-1.1957542) q[1];
sx q[1];
rz(3.0072838) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90420049) q[0];
sx q[0];
rz(-0.89590329) q[0];
sx q[0];
rz(-0.86752059) q[0];
x q[1];
rz(-2.9555126) q[2];
sx q[2];
rz(-1.1717516) q[2];
sx q[2];
rz(-3.0732791) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1301535) q[1];
sx q[1];
rz(-1.7053889) q[1];
sx q[1];
rz(-2.9573739) q[1];
rz(-0.28388309) q[3];
sx q[3];
rz(-1.1177269) q[3];
sx q[3];
rz(2.4456152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.40992752) q[2];
sx q[2];
rz(-1.8525367) q[2];
sx q[2];
rz(-0.64201391) q[2];
rz(2.0783453) q[3];
sx q[3];
rz(-0.65170538) q[3];
sx q[3];
rz(-1.0412019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6006271) q[0];
sx q[0];
rz(-0.0890812) q[0];
sx q[0];
rz(0.11216057) q[0];
rz(1.3345831) q[1];
sx q[1];
rz(-0.59252512) q[1];
sx q[1];
rz(2.1122011) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6717259) q[0];
sx q[0];
rz(-0.70505667) q[0];
sx q[0];
rz(1.6204349) q[0];
x q[1];
rz(1.5310982) q[2];
sx q[2];
rz(-1.8744812) q[2];
sx q[2];
rz(-0.025321753) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.411602) q[1];
sx q[1];
rz(-1.3407602) q[1];
sx q[1];
rz(0.045129808) q[1];
rz(-pi) q[2];
rz(-2.6736027) q[3];
sx q[3];
rz(-1.6525998) q[3];
sx q[3];
rz(1.8379267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7182497) q[2];
sx q[2];
rz(-0.074904718) q[2];
sx q[2];
rz(-3.1035799) q[2];
rz(0.612261) q[3];
sx q[3];
rz(-2.2050048) q[3];
sx q[3];
rz(3.0003701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24899471) q[0];
sx q[0];
rz(-1.4709512) q[0];
sx q[0];
rz(0.41879642) q[0];
rz(1.2576125) q[1];
sx q[1];
rz(-1.6323171) q[1];
sx q[1];
rz(2.9990101) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.044238) q[0];
sx q[0];
rz(-1.729106) q[0];
sx q[0];
rz(1.5990785) q[0];
rz(1.293574) q[2];
sx q[2];
rz(-0.96053329) q[2];
sx q[2];
rz(-1.1278314) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.31227641) q[1];
sx q[1];
rz(-1.3061151) q[1];
sx q[1];
rz(-2.1618202) q[1];
rz(-pi) q[2];
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
rz(1.3451781) q[2];
sx q[2];
rz(-0.08064457) q[2];
sx q[2];
rz(-0.26819116) q[2];
rz(1.7372519) q[3];
sx q[3];
rz(-2.0495448) q[3];
sx q[3];
rz(-2.0973189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0161229) q[0];
sx q[0];
rz(-0.37624993) q[0];
sx q[0];
rz(-2.7952623) q[0];
rz(-0.12915962) q[1];
sx q[1];
rz(-1.8864514) q[1];
sx q[1];
rz(1.4056978) q[1];
rz(0.62025537) q[2];
sx q[2];
rz(-0.56369416) q[2];
sx q[2];
rz(-0.77947215) q[2];
rz(-0.066349647) q[3];
sx q[3];
rz(-1.8649615) q[3];
sx q[3];
rz(-2.1129114) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
