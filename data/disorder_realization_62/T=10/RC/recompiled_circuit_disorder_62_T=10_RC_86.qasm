OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9010889) q[0];
sx q[0];
rz(-0.17740372) q[0];
sx q[0];
rz(1.1343962) q[0];
rz(-1.9534684) q[1];
sx q[1];
rz(-1.0367353) q[1];
sx q[1];
rz(-2.477975) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9288869) q[0];
sx q[0];
rz(-2.0318673) q[0];
sx q[0];
rz(-1.5443718) q[0];
rz(2.7825836) q[2];
sx q[2];
rz(-0.81072545) q[2];
sx q[2];
rz(0.63149482) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.692688) q[1];
sx q[1];
rz(-1.2346039) q[1];
sx q[1];
rz(-0.2775788) q[1];
rz(-1.1575559) q[3];
sx q[3];
rz(-2.8765656) q[3];
sx q[3];
rz(2.7812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2628281) q[2];
sx q[2];
rz(-0.44439134) q[2];
sx q[2];
rz(0.051068548) q[2];
rz(0.55705327) q[3];
sx q[3];
rz(-0.80018187) q[3];
sx q[3];
rz(-1.5867656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58650815) q[0];
sx q[0];
rz(-0.83207911) q[0];
sx q[0];
rz(2.5449975) q[0];
rz(2.3157628) q[1];
sx q[1];
rz(-1.700371) q[1];
sx q[1];
rz(-1.9155496) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37218371) q[0];
sx q[0];
rz(-0.49159494) q[0];
sx q[0];
rz(-2.7093637) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77483564) q[2];
sx q[2];
rz(-0.92445395) q[2];
sx q[2];
rz(-1.6006084) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5920168) q[1];
sx q[1];
rz(-1.5703652) q[1];
sx q[1];
rz(1.3509343) q[1];
x q[2];
rz(3.1321208) q[3];
sx q[3];
rz(-0.99820271) q[3];
sx q[3];
rz(3.115311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3423959) q[2];
sx q[2];
rz(-1.1725972) q[2];
sx q[2];
rz(-0.33102316) q[2];
rz(-2.3349169) q[3];
sx q[3];
rz(-1.2253864) q[3];
sx q[3];
rz(-1.7139009) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1598635) q[0];
sx q[0];
rz(-1.230343) q[0];
sx q[0];
rz(1.8925517) q[0];
rz(3.0535835) q[1];
sx q[1];
rz(-2.0188589) q[1];
sx q[1];
rz(2.1121315) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91988504) q[0];
sx q[0];
rz(-1.789933) q[0];
sx q[0];
rz(2.1973781) q[0];
rz(1.8573895) q[2];
sx q[2];
rz(-0.9409875) q[2];
sx q[2];
rz(2.0558002) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6368235) q[1];
sx q[1];
rz(-1.4154589) q[1];
sx q[1];
rz(1.7117281) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73996468) q[3];
sx q[3];
rz(-1.660534) q[3];
sx q[3];
rz(1.2147853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.007894667) q[2];
sx q[2];
rz(-1.7241314) q[2];
sx q[2];
rz(-0.50764817) q[2];
rz(-1.7525904) q[3];
sx q[3];
rz(-0.32998431) q[3];
sx q[3];
rz(1.9784137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4841109) q[0];
sx q[0];
rz(-2.8634475) q[0];
sx q[0];
rz(1.5456276) q[0];
rz(-1.0428628) q[1];
sx q[1];
rz(-1.9680126) q[1];
sx q[1];
rz(1.625659) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0632616) q[0];
sx q[0];
rz(-2.445979) q[0];
sx q[0];
rz(-0.22380933) q[0];
x q[1];
rz(-1.707294) q[2];
sx q[2];
rz(-0.56328661) q[2];
sx q[2];
rz(2.0493281) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8856636) q[1];
sx q[1];
rz(-0.71223488) q[1];
sx q[1];
rz(-2.7182012) q[1];
x q[2];
rz(-1.3514148) q[3];
sx q[3];
rz(-0.99321584) q[3];
sx q[3];
rz(2.1882309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3778014) q[2];
sx q[2];
rz(-0.8447454) q[2];
sx q[2];
rz(1.2949004) q[2];
rz(-3.0002248) q[3];
sx q[3];
rz(-2.5989792) q[3];
sx q[3];
rz(-2.1550089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35448733) q[0];
sx q[0];
rz(-1.1802477) q[0];
sx q[0];
rz(-1.6217344) q[0];
rz(-0.63201085) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(0.89486665) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7475815) q[0];
sx q[0];
rz(-2.3590439) q[0];
sx q[0];
rz(0.6015908) q[0];
rz(-pi) q[1];
rz(-0.72563719) q[2];
sx q[2];
rz(-1.1302395) q[2];
sx q[2];
rz(2.5051136) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.27069651) q[1];
sx q[1];
rz(-2.0082698) q[1];
sx q[1];
rz(2.2938964) q[1];
rz(0.026167913) q[3];
sx q[3];
rz(-2.8550365) q[3];
sx q[3];
rz(-2.9122796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.616509) q[2];
sx q[2];
rz(-0.36642763) q[2];
sx q[2];
rz(-1.9990702) q[2];
rz(-0.74674314) q[3];
sx q[3];
rz(-1.8871504) q[3];
sx q[3];
rz(-2.0660627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58364761) q[0];
sx q[0];
rz(-0.32145158) q[0];
sx q[0];
rz(1.7640132) q[0];
rz(-0.47239834) q[1];
sx q[1];
rz(-0.51858416) q[1];
sx q[1];
rz(-2.6766052) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2705921) q[0];
sx q[0];
rz(-2.5809079) q[0];
sx q[0];
rz(-2.2339348) q[0];
rz(-pi) q[1];
rz(-0.34157413) q[2];
sx q[2];
rz(-1.5014868) q[2];
sx q[2];
rz(0.14788936) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2496693) q[1];
sx q[1];
rz(-1.6185456) q[1];
sx q[1];
rz(0.93530099) q[1];
rz(-1.7808077) q[3];
sx q[3];
rz(-0.87267733) q[3];
sx q[3];
rz(-0.39892808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.650699) q[2];
sx q[2];
rz(-1.2630254) q[2];
sx q[2];
rz(-0.4450376) q[2];
rz(-2.2079091) q[3];
sx q[3];
rz(-1.7088339) q[3];
sx q[3];
rz(-2.8745108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12748195) q[0];
sx q[0];
rz(-1.575379) q[0];
sx q[0];
rz(1.7215464) q[0];
rz(-0.02380112) q[1];
sx q[1];
rz(-2.5271466) q[1];
sx q[1];
rz(2.9856317) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1056846) q[0];
sx q[0];
rz(-1.3996291) q[0];
sx q[0];
rz(-0.25733421) q[0];
x q[1];
rz(2.4620352) q[2];
sx q[2];
rz(-0.62629269) q[2];
sx q[2];
rz(1.2223513) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.012688795) q[1];
sx q[1];
rz(-1.745599) q[1];
sx q[1];
rz(-0.054794475) q[1];
rz(-2.0701253) q[3];
sx q[3];
rz(-1.0896177) q[3];
sx q[3];
rz(1.037998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0722787) q[2];
sx q[2];
rz(-2.5645655) q[2];
sx q[2];
rz(2.0689266) q[2];
rz(2.8159451) q[3];
sx q[3];
rz(-1.1762534) q[3];
sx q[3];
rz(1.6252888) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1919365) q[0];
sx q[0];
rz(-2.7392116) q[0];
sx q[0];
rz(-2.8038213) q[0];
rz(-1.0900963) q[1];
sx q[1];
rz(-0.97507674) q[1];
sx q[1];
rz(-0.24857323) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08911207) q[0];
sx q[0];
rz(-0.57894527) q[0];
sx q[0];
rz(0.46778932) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0806662) q[2];
sx q[2];
rz(-1.14398) q[2];
sx q[2];
rz(-0.61444297) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.92330248) q[1];
sx q[1];
rz(-2.3308838) q[1];
sx q[1];
rz(0.87358012) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.076898889) q[3];
sx q[3];
rz(-1.4014763) q[3];
sx q[3];
rz(-2.3016735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8481855) q[2];
sx q[2];
rz(-2.6082787) q[2];
sx q[2];
rz(3.0548813) q[2];
rz(-0.48197204) q[3];
sx q[3];
rz(-0.50589365) q[3];
sx q[3];
rz(1.1530676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50865737) q[0];
sx q[0];
rz(-2.1338699) q[0];
sx q[0];
rz(-0.28276643) q[0];
rz(2.4400318) q[1];
sx q[1];
rz(-2.3200254) q[1];
sx q[1];
rz(1.823002) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.628172) q[0];
sx q[0];
rz(-1.6512198) q[0];
sx q[0];
rz(2.0128987) q[0];
x q[1];
rz(0.65638541) q[2];
sx q[2];
rz(-2.0771386) q[2];
sx q[2];
rz(0.93751794) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1773771) q[1];
sx q[1];
rz(-2.4417158) q[1];
sx q[1];
rz(1.3436951) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3913279) q[3];
sx q[3];
rz(-1.2816458) q[3];
sx q[3];
rz(-0.17444785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41137722) q[2];
sx q[2];
rz(-0.29835478) q[2];
sx q[2];
rz(-0.39917699) q[2];
rz(2.2579851) q[3];
sx q[3];
rz(-1.4727605) q[3];
sx q[3];
rz(1.9201027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3783962) q[0];
sx q[0];
rz(-2.3662687) q[0];
sx q[0];
rz(1.5426853) q[0];
rz(2.0762766) q[1];
sx q[1];
rz(-0.97384802) q[1];
sx q[1];
rz(-1.261196) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48860088) q[0];
sx q[0];
rz(-1.6761259) q[0];
sx q[0];
rz(-0.41098849) q[0];
rz(-2.5433259) q[2];
sx q[2];
rz(-1.5891799) q[2];
sx q[2];
rz(0.53668864) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.8775455) q[1];
sx q[1];
rz(-2.0685158) q[1];
sx q[1];
rz(-1.4977786) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.946978) q[3];
sx q[3];
rz(-2.448304) q[3];
sx q[3];
rz(0.64893901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5836872) q[2];
sx q[2];
rz(-0.62290278) q[2];
sx q[2];
rz(-2.5718001) q[2];
rz(1.9231046) q[3];
sx q[3];
rz(-2.4945538) q[3];
sx q[3];
rz(-2.5789554) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8284843) q[0];
sx q[0];
rz(-0.88043558) q[0];
sx q[0];
rz(1.2783929) q[0];
rz(-0.60824153) q[1];
sx q[1];
rz(-0.47641644) q[1];
sx q[1];
rz(0.48412916) q[1];
rz(3.1153395) q[2];
sx q[2];
rz(-0.99925169) q[2];
sx q[2];
rz(3.0796438) q[2];
rz(0.59109296) q[3];
sx q[3];
rz(-1.686284) q[3];
sx q[3];
rz(-1.4808663) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
