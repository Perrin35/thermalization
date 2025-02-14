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
rz(-0.64088351) q[1];
sx q[1];
rz(-1.0008608) q[1];
sx q[1];
rz(0.34520087) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76287718) q[0];
sx q[0];
rz(-1.595531) q[0];
sx q[0];
rz(-1.6616761) q[0];
rz(-pi) q[1];
rz(2.7852374) q[2];
sx q[2];
rz(-0.40268597) q[2];
sx q[2];
rz(2.8860983) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4507323) q[1];
sx q[1];
rz(-1.9885411) q[1];
sx q[1];
rz(-0.91207204) q[1];
x q[2];
rz(-3.1253417) q[3];
sx q[3];
rz(-1.645184) q[3];
sx q[3];
rz(2.1824238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9125646) q[2];
sx q[2];
rz(-1.8077069) q[2];
sx q[2];
rz(0.20797569) q[2];
rz(0.29911706) q[3];
sx q[3];
rz(-0.57204539) q[3];
sx q[3];
rz(2.0007029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3308554) q[0];
sx q[0];
rz(-2.9006697) q[0];
sx q[0];
rz(-0.99288565) q[0];
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
rz(-2.3925288) q[0];
sx q[0];
rz(-1.3925902) q[0];
sx q[0];
rz(-3.130413) q[0];
rz(-pi) q[1];
rz(1.0335017) q[2];
sx q[2];
rz(-0.8391434) q[2];
sx q[2];
rz(-3.0947859) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.44821139) q[1];
sx q[1];
rz(-0.44084545) q[1];
sx q[1];
rz(-0.54852672) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1466188) q[3];
sx q[3];
rz(-2.6335003) q[3];
sx q[3];
rz(1.0095694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2449067) q[2];
sx q[2];
rz(-1.2865571) q[2];
sx q[2];
rz(-2.1265105) q[2];
rz(-0.24909881) q[3];
sx q[3];
rz(-2.2738012) q[3];
sx q[3];
rz(-2.7740313) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2772813) q[0];
sx q[0];
rz(-0.69121498) q[0];
sx q[0];
rz(-0.15750289) q[0];
rz(2.1306254) q[1];
sx q[1];
rz(-2.8087661) q[1];
sx q[1];
rz(-0.13883042) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.39108) q[0];
sx q[0];
rz(-1.2475999) q[0];
sx q[0];
rz(1.9531519) q[0];
rz(1.9857668) q[2];
sx q[2];
rz(-0.89902821) q[2];
sx q[2];
rz(-0.49723724) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0142747) q[1];
sx q[1];
rz(-1.7671314) q[1];
sx q[1];
rz(-1.8683968) q[1];
rz(2.8582358) q[3];
sx q[3];
rz(-1.4025926) q[3];
sx q[3];
rz(3.0137872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6802754) q[2];
sx q[2];
rz(-2.2291144) q[2];
sx q[2];
rz(2.7834564) q[2];
rz(-0.67874587) q[3];
sx q[3];
rz(-0.94921422) q[3];
sx q[3];
rz(2.1024735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.045227483) q[0];
sx q[0];
rz(-0.97310936) q[0];
sx q[0];
rz(2.8367693) q[0];
rz(0.40799704) q[1];
sx q[1];
rz(-1.4381189) q[1];
sx q[1];
rz(-2.2136484) q[1];
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
x q[1];
rz(-0.94522743) q[2];
sx q[2];
rz(-2.0778227) q[2];
sx q[2];
rz(-0.27238174) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2603307) q[1];
sx q[1];
rz(-1.7219543) q[1];
sx q[1];
rz(2.4511454) q[1];
x q[2];
rz(-2.4797012) q[3];
sx q[3];
rz(-2.3958979) q[3];
sx q[3];
rz(-2.6357366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1912332) q[2];
sx q[2];
rz(-0.57098907) q[2];
sx q[2];
rz(-2.9534269) q[2];
rz(1.2763216) q[3];
sx q[3];
rz(-1.8264344) q[3];
sx q[3];
rz(-1.4307384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6998049) q[0];
sx q[0];
rz(-2.7975174) q[0];
sx q[0];
rz(3.1222043) q[0];
rz(1.060932) q[1];
sx q[1];
rz(-1.414199) q[1];
sx q[1];
rz(-1.9020938) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083154924) q[0];
sx q[0];
rz(-2.4920336) q[0];
sx q[0];
rz(-2.044528) q[0];
x q[1];
rz(2.4475736) q[2];
sx q[2];
rz(-0.87022129) q[2];
sx q[2];
rz(-1.7050336) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.13889276) q[1];
sx q[1];
rz(-2.5475916) q[1];
sx q[1];
rz(0.51704375) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.89770384) q[3];
sx q[3];
rz(-0.60808676) q[3];
sx q[3];
rz(1.2230108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0197319) q[2];
sx q[2];
rz(-1.8283565) q[2];
sx q[2];
rz(-1.8149553) q[2];
rz(-0.2462247) q[3];
sx q[3];
rz(-2.0387869) q[3];
sx q[3];
rz(-2.2897913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5354079) q[0];
sx q[0];
rz(-0.98537213) q[0];
sx q[0];
rz(-0.73915172) q[0];
rz(0.34573653) q[1];
sx q[1];
rz(-1.3290936) q[1];
sx q[1];
rz(-0.87127042) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5369878) q[0];
sx q[0];
rz(-0.78442998) q[0];
sx q[0];
rz(-3.0434199) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61001444) q[2];
sx q[2];
rz(-2.2918309) q[2];
sx q[2];
rz(0.80243669) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4193486) q[1];
sx q[1];
rz(-2.7461395) q[1];
sx q[1];
rz(-1.761318) q[1];
rz(0.61136758) q[3];
sx q[3];
rz(-1.1725559) q[3];
sx q[3];
rz(0.37722019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.31200108) q[2];
sx q[2];
rz(-2.5025554) q[2];
sx q[2];
rz(-0.87257067) q[2];
rz(1.568659) q[3];
sx q[3];
rz(-0.60397732) q[3];
sx q[3];
rz(2.2085371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0912112) q[0];
sx q[0];
rz(-2.8822883) q[0];
sx q[0];
rz(-0.76164371) q[0];
rz(-0.42539445) q[1];
sx q[1];
rz(-1.1401221) q[1];
sx q[1];
rz(-0.92686191) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2970718) q[0];
sx q[0];
rz(-1.8064587) q[0];
sx q[0];
rz(-1.8051487) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4285467) q[2];
sx q[2];
rz(-2.2879061) q[2];
sx q[2];
rz(2.6756659) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1897404) q[1];
sx q[1];
rz(-1.2092672) q[1];
sx q[1];
rz(0.23747634) q[1];
rz(-pi) q[2];
x q[2];
rz(0.94759946) q[3];
sx q[3];
rz(-1.609908) q[3];
sx q[3];
rz(-0.87621237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0285792) q[2];
sx q[2];
rz(-0.69602746) q[2];
sx q[2];
rz(-2.2704303) q[2];
rz(-0.29843676) q[3];
sx q[3];
rz(-1.2454183) q[3];
sx q[3];
rz(-1.3309853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4153862) q[0];
sx q[0];
rz(-2.6415249) q[0];
sx q[0];
rz(-2.723208) q[0];
rz(0.2977953) q[1];
sx q[1];
rz(-1.1957542) q[1];
sx q[1];
rz(0.13430886) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1100562) q[0];
sx q[0];
rz(-2.2087065) q[0];
sx q[0];
rz(0.67968926) q[0];
rz(-1.9842582) q[2];
sx q[2];
rz(-2.7034139) q[2];
sx q[2];
rz(0.38288051) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.93523504) q[1];
sx q[1];
rz(-2.9138953) q[1];
sx q[1];
rz(2.5050194) q[1];
rz(-pi) q[2];
rz(-2.8577096) q[3];
sx q[3];
rz(-2.0238658) q[3];
sx q[3];
rz(2.4456152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.40992752) q[2];
sx q[2];
rz(-1.2890559) q[2];
sx q[2];
rz(-0.64201391) q[2];
rz(2.0783453) q[3];
sx q[3];
rz(-2.4898873) q[3];
sx q[3];
rz(1.0412019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5409656) q[0];
sx q[0];
rz(-3.0525115) q[0];
sx q[0];
rz(-3.0294321) q[0];
rz(-1.3345831) q[1];
sx q[1];
rz(-2.5490675) q[1];
sx q[1];
rz(2.1122011) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4698668) q[0];
sx q[0];
rz(-2.436536) q[0];
sx q[0];
rz(-1.6204349) q[0];
x q[1];
rz(-0.12597398) q[2];
sx q[2];
rz(-2.8354037) q[2];
sx q[2];
rz(0.15737113) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.411602) q[1];
sx q[1];
rz(-1.3407602) q[1];
sx q[1];
rz(-0.045129808) q[1];
rz(0.46799) q[3];
sx q[3];
rz(-1.6525998) q[3];
sx q[3];
rz(1.8379267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7182497) q[2];
sx q[2];
rz(-3.0666879) q[2];
sx q[2];
rz(-3.1035799) q[2];
rz(0.612261) q[3];
sx q[3];
rz(-2.2050048) q[3];
sx q[3];
rz(-0.14122252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24899471) q[0];
sx q[0];
rz(-1.6706415) q[0];
sx q[0];
rz(0.41879642) q[0];
rz(-1.8839802) q[1];
sx q[1];
rz(-1.5092756) q[1];
sx q[1];
rz(0.14258252) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.044238) q[0];
sx q[0];
rz(-1.729106) q[0];
sx q[0];
rz(1.5425141) q[0];
rz(-pi) q[1];
rz(2.5129287) q[2];
sx q[2];
rz(-1.7970048) q[2];
sx q[2];
rz(-2.5369801) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8293162) q[1];
sx q[1];
rz(-1.8354776) q[1];
sx q[1];
rz(-2.1618202) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0311277) q[3];
sx q[3];
rz(-1.2602196) q[3];
sx q[3];
rz(-2.6643857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3451781) q[2];
sx q[2];
rz(-0.08064457) q[2];
sx q[2];
rz(-2.8734015) q[2];
rz(1.7372519) q[3];
sx q[3];
rz(-1.0920478) q[3];
sx q[3];
rz(-1.0442737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1254697) q[0];
sx q[0];
rz(-0.37624993) q[0];
sx q[0];
rz(-2.7952623) q[0];
rz(-0.12915962) q[1];
sx q[1];
rz(-1.8864514) q[1];
sx q[1];
rz(1.4056978) q[1];
rz(2.5213373) q[2];
sx q[2];
rz(-2.5778985) q[2];
sx q[2];
rz(2.3621205) q[2];
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
