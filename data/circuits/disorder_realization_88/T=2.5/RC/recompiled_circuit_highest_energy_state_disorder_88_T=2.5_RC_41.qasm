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
rz(-3.1336194) q[0];
sx q[0];
rz(-2.70533) q[0];
sx q[0];
rz(-2.1313957) q[0];
rz(0.17207347) q[1];
sx q[1];
rz(-2.3354524) q[1];
sx q[1];
rz(-1.0239209) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9925594) q[0];
sx q[0];
rz(-1.5776538) q[0];
sx q[0];
rz(-1.0829145) q[0];
x q[1];
rz(-1.3390117) q[2];
sx q[2];
rz(-1.9531774) q[2];
sx q[2];
rz(0.2222375) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.55320239) q[1];
sx q[1];
rz(-1.4103048) q[1];
sx q[1];
rz(-3.1218916) q[1];
rz(-pi) q[2];
rz(0.2497345) q[3];
sx q[3];
rz(-1.485255) q[3];
sx q[3];
rz(-1.5126603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.94850928) q[2];
sx q[2];
rz(-0.32863363) q[2];
sx q[2];
rz(-2.9060717) q[2];
rz(1.0282907) q[3];
sx q[3];
rz(-1.2582658) q[3];
sx q[3];
rz(-1.9839015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.704772) q[0];
sx q[0];
rz(-1.3587767) q[0];
sx q[0];
rz(-0.92192465) q[0];
rz(3.0251265) q[1];
sx q[1];
rz(-1.4215697) q[1];
sx q[1];
rz(2.2540653) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5619316) q[0];
sx q[0];
rz(-1.4401739) q[0];
sx q[0];
rz(-1.3829253) q[0];
rz(-2.9824184) q[2];
sx q[2];
rz(-2.7870057) q[2];
sx q[2];
rz(1.872242) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0640338) q[1];
sx q[1];
rz(-2.56018) q[1];
sx q[1];
rz(-0.20661776) q[1];
rz(-pi) q[2];
rz(1.4027897) q[3];
sx q[3];
rz(-0.6989494) q[3];
sx q[3];
rz(-0.93937554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.55887115) q[2];
sx q[2];
rz(-1.0251986) q[2];
sx q[2];
rz(0.23898807) q[2];
rz(-0.72238266) q[3];
sx q[3];
rz(-2.3014849) q[3];
sx q[3];
rz(1.1649789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8257985) q[0];
sx q[0];
rz(-2.1738985) q[0];
sx q[0];
rz(-2.353299) q[0];
rz(-1.7514508) q[1];
sx q[1];
rz(-2.7056521) q[1];
sx q[1];
rz(1.4424666) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6955439) q[0];
sx q[0];
rz(-0.15803495) q[0];
sx q[0];
rz(-2.0331073) q[0];
rz(-pi) q[1];
x q[1];
rz(1.945247) q[2];
sx q[2];
rz(-1.8467287) q[2];
sx q[2];
rz(0.79513136) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4462112) q[1];
sx q[1];
rz(-1.3331116) q[1];
sx q[1];
rz(2.2878245) q[1];
rz(-pi) q[2];
rz(0.17876153) q[3];
sx q[3];
rz(-1.074665) q[3];
sx q[3];
rz(2.2876491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73551377) q[2];
sx q[2];
rz(-1.9903851) q[2];
sx q[2];
rz(2.8847412) q[2];
rz(2.2779321) q[3];
sx q[3];
rz(-2.05859) q[3];
sx q[3];
rz(-2.5679722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0470806) q[0];
sx q[0];
rz(-0.032945078) q[0];
sx q[0];
rz(1.5091913) q[0];
rz(2.8250561) q[1];
sx q[1];
rz(-1.5959975) q[1];
sx q[1];
rz(-2.1260156) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2779396) q[0];
sx q[0];
rz(-2.4653593) q[0];
sx q[0];
rz(0.77267026) q[0];
x q[1];
rz(-0.79358856) q[2];
sx q[2];
rz(-1.3745135) q[2];
sx q[2];
rz(-1.0283141) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9719661) q[1];
sx q[1];
rz(-1.7952654) q[1];
sx q[1];
rz(-1.1826993) q[1];
rz(1.7523132) q[3];
sx q[3];
rz(-2.3861775) q[3];
sx q[3];
rz(0.28873539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.29869002) q[2];
sx q[2];
rz(-2.0899253) q[2];
sx q[2];
rz(-1.6896557) q[2];
rz(-1.2339833) q[3];
sx q[3];
rz(-2.0472287) q[3];
sx q[3];
rz(0.88172495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3877617) q[0];
sx q[0];
rz(-1.8330638) q[0];
sx q[0];
rz(2.0552788) q[0];
rz(-1.7408675) q[1];
sx q[1];
rz(-0.81175214) q[1];
sx q[1];
rz(0.0033671826) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8400492) q[0];
sx q[0];
rz(-2.1652984) q[0];
sx q[0];
rz(-1.5778753) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0900159) q[2];
sx q[2];
rz(-0.64489105) q[2];
sx q[2];
rz(0.99586785) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0778179) q[1];
sx q[1];
rz(-1.5146202) q[1];
sx q[1];
rz(0.50775524) q[1];
rz(-pi) q[2];
rz(-1.5448872) q[3];
sx q[3];
rz(-0.97086519) q[3];
sx q[3];
rz(2.5480218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9747666) q[2];
sx q[2];
rz(-2.4415015) q[2];
sx q[2];
rz(-0.33494803) q[2];
rz(2.8454928) q[3];
sx q[3];
rz(-0.94303232) q[3];
sx q[3];
rz(-0.12224841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5122546) q[0];
sx q[0];
rz(-0.4991931) q[0];
sx q[0];
rz(-2.7235624) q[0];
rz(2.4755075) q[1];
sx q[1];
rz(-1.2434375) q[1];
sx q[1];
rz(-0.24806771) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92221602) q[0];
sx q[0];
rz(-1.9016978) q[0];
sx q[0];
rz(-1.0537345) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.054912189) q[2];
sx q[2];
rz(-1.1217199) q[2];
sx q[2];
rz(2.1841846) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.848143) q[1];
sx q[1];
rz(-2.0781233) q[1];
sx q[1];
rz(0.44815843) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8855528) q[3];
sx q[3];
rz(-2.1845792) q[3];
sx q[3];
rz(-2.7308794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.69998133) q[2];
sx q[2];
rz(-1.8786083) q[2];
sx q[2];
rz(1.9786394) q[2];
rz(0.59398389) q[3];
sx q[3];
rz(-1.5409527) q[3];
sx q[3];
rz(-2.1514413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6043337) q[0];
sx q[0];
rz(-2.8626677) q[0];
sx q[0];
rz(3.0766686) q[0];
rz(0.36554947) q[1];
sx q[1];
rz(-1.7100916) q[1];
sx q[1];
rz(-2.4117267) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6014272) q[0];
sx q[0];
rz(-1.6807204) q[0];
sx q[0];
rz(2.2137627) q[0];
rz(-1.3781629) q[2];
sx q[2];
rz(-0.81708497) q[2];
sx q[2];
rz(-1.1258923) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5090885) q[1];
sx q[1];
rz(-1.4756241) q[1];
sx q[1];
rz(-1.0126963) q[1];
rz(-2.4555931) q[3];
sx q[3];
rz(-1.2427875) q[3];
sx q[3];
rz(2.6552116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2713752) q[2];
sx q[2];
rz(-1.368618) q[2];
sx q[2];
rz(1.082487) q[2];
rz(1.6535951) q[3];
sx q[3];
rz(-0.36181417) q[3];
sx q[3];
rz(-0.35058072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7545886) q[0];
sx q[0];
rz(-0.8828187) q[0];
sx q[0];
rz(0.3311232) q[0];
rz(1.4777199) q[1];
sx q[1];
rz(-0.96550566) q[1];
sx q[1];
rz(-0.32041916) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5950111) q[0];
sx q[0];
rz(-0.10215952) q[0];
sx q[0];
rz(-0.25545303) q[0];
x q[1];
rz(1.4244979) q[2];
sx q[2];
rz(-1.723389) q[2];
sx q[2];
rz(2.1426147) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.85727188) q[1];
sx q[1];
rz(-1.1445938) q[1];
sx q[1];
rz(1.5537474) q[1];
x q[2];
rz(-1.6410758) q[3];
sx q[3];
rz(-0.68551862) q[3];
sx q[3];
rz(0.34834114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.42935818) q[2];
sx q[2];
rz(-1.7244491) q[2];
sx q[2];
rz(-2.7799535) q[2];
rz(-2.2278191) q[3];
sx q[3];
rz(-2.8445966) q[3];
sx q[3];
rz(-0.44338068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(-0.2509895) q[0];
sx q[0];
rz(-0.78762233) q[0];
sx q[0];
rz(0.11601624) q[0];
rz(1.3277671) q[1];
sx q[1];
rz(-1.1632183) q[1];
sx q[1];
rz(-1.2001002) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3108771) q[0];
sx q[0];
rz(-2.9492464) q[0];
sx q[0];
rz(-0.098523454) q[0];
rz(0.88516) q[2];
sx q[2];
rz(-1.4421534) q[2];
sx q[2];
rz(-2.7264915) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1397651) q[1];
sx q[1];
rz(-1.5854302) q[1];
sx q[1];
rz(-2.57994) q[1];
rz(-pi) q[2];
rz(0.16606826) q[3];
sx q[3];
rz(-1.8072053) q[3];
sx q[3];
rz(1.1354367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8414595) q[2];
sx q[2];
rz(-1.7068784) q[2];
sx q[2];
rz(2.8099698) q[2];
rz(-0.2837818) q[3];
sx q[3];
rz(-2.0399358) q[3];
sx q[3];
rz(2.7200429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3109741) q[0];
sx q[0];
rz(-2.0299439) q[0];
sx q[0];
rz(-2.9750138) q[0];
rz(-0.84959787) q[1];
sx q[1];
rz(-2.6764937) q[1];
sx q[1];
rz(-0.073624484) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66921418) q[0];
sx q[0];
rz(-1.3166691) q[0];
sx q[0];
rz(0.73331244) q[0];
rz(2.6643986) q[2];
sx q[2];
rz(-1.0721803) q[2];
sx q[2];
rz(-2.3248364) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.4576118) q[1];
sx q[1];
rz(-1.906412) q[1];
sx q[1];
rz(2.6102351) q[1];
rz(-2.5166148) q[3];
sx q[3];
rz(-2.4132937) q[3];
sx q[3];
rz(2.7788904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.9061884) q[2];
sx q[2];
rz(-0.81934682) q[2];
sx q[2];
rz(-0.6582312) q[2];
rz(1.341358) q[3];
sx q[3];
rz(-0.96786371) q[3];
sx q[3];
rz(-0.85097504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.3157208) q[0];
sx q[0];
rz(-0.82875874) q[0];
sx q[0];
rz(2.2525633) q[0];
rz(2.6201164) q[1];
sx q[1];
rz(-1.5670525) q[1];
sx q[1];
rz(-1.570931) q[1];
rz(-2.1076219) q[2];
sx q[2];
rz(-2.2848717) q[2];
sx q[2];
rz(0.90938766) q[2];
rz(0.23052774) q[3];
sx q[3];
rz(-2.5067825) q[3];
sx q[3];
rz(1.4092177) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
