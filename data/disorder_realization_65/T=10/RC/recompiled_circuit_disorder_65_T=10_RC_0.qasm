OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31057519) q[0];
sx q[0];
rz(-1.0520881) q[0];
sx q[0];
rz(-1.4927827) q[0];
rz(0.86039034) q[1];
sx q[1];
rz(-2.4465423) q[1];
sx q[1];
rz(2.066943) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63700919) q[0];
sx q[0];
rz(-1.552581) q[0];
sx q[0];
rz(-1.5760742) q[0];
rz(-pi) q[1];
rz(-1.4104112) q[2];
sx q[2];
rz(-0.78750247) q[2];
sx q[2];
rz(-3.0410142) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5584471) q[1];
sx q[1];
rz(-1.4902643) q[1];
sx q[1];
rz(0.12195485) q[1];
rz(-0.53332897) q[3];
sx q[3];
rz(-1.9682069) q[3];
sx q[3];
rz(-1.0226137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3124714) q[2];
sx q[2];
rz(-2.1437936) q[2];
sx q[2];
rz(-1.5134229) q[2];
rz(1.9681905) q[3];
sx q[3];
rz(-1.4459926) q[3];
sx q[3];
rz(0.2127969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5728977) q[0];
sx q[0];
rz(-0.64210367) q[0];
sx q[0];
rz(-1.9182385) q[0];
rz(2.9808295) q[1];
sx q[1];
rz(-1.5784966) q[1];
sx q[1];
rz(1.5011903) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0813635) q[0];
sx q[0];
rz(-1.4481359) q[0];
sx q[0];
rz(1.5958022) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8295849) q[2];
sx q[2];
rz(-1.1216251) q[2];
sx q[2];
rz(-2.1829407) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7306819) q[1];
sx q[1];
rz(-0.69697471) q[1];
sx q[1];
rz(-0.26941401) q[1];
rz(3.1355751) q[3];
sx q[3];
rz(-1.179751) q[3];
sx q[3];
rz(-1.0139549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.286065) q[2];
sx q[2];
rz(-2.3748886) q[2];
sx q[2];
rz(-1.6274874) q[2];
rz(1.9747915) q[3];
sx q[3];
rz(-1.6191926) q[3];
sx q[3];
rz(-0.91439009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0062155) q[0];
sx q[0];
rz(-1.051798) q[0];
sx q[0];
rz(-1.0789543) q[0];
rz(-1.9619933) q[1];
sx q[1];
rz(-1.0410573) q[1];
sx q[1];
rz(-1.5037781) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30341879) q[0];
sx q[0];
rz(-1.5056599) q[0];
sx q[0];
rz(-1.6584048) q[0];
rz(2.1554699) q[2];
sx q[2];
rz(-1.1650656) q[2];
sx q[2];
rz(-3.0509146) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1551731) q[1];
sx q[1];
rz(-1.5281614) q[1];
sx q[1];
rz(-2.8405872) q[1];
rz(-pi) q[2];
rz(1.740326) q[3];
sx q[3];
rz(-2.4781514) q[3];
sx q[3];
rz(-1.1273813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.43061313) q[2];
sx q[2];
rz(-0.47982275) q[2];
sx q[2];
rz(0.7437931) q[2];
rz(-0.25742325) q[3];
sx q[3];
rz(-1.0835203) q[3];
sx q[3];
rz(-2.553489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1733615) q[0];
sx q[0];
rz(-3.0259202) q[0];
sx q[0];
rz(-2.0898576) q[0];
rz(-2.5412718) q[1];
sx q[1];
rz(-0.61906639) q[1];
sx q[1];
rz(1.1490885) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7950492) q[0];
sx q[0];
rz(-1.8115037) q[0];
sx q[0];
rz(1.6853149) q[0];
x q[1];
rz(-0.43196584) q[2];
sx q[2];
rz(-1.0216121) q[2];
sx q[2];
rz(-0.84749046) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9149575) q[1];
sx q[1];
rz(-1.8640222) q[1];
sx q[1];
rz(0.83421591) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37255128) q[3];
sx q[3];
rz(-0.77790341) q[3];
sx q[3];
rz(-2.3140964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8950618) q[2];
sx q[2];
rz(-1.4738513) q[2];
sx q[2];
rz(-2.4728298) q[2];
rz(1.8956005) q[3];
sx q[3];
rz(-0.99530882) q[3];
sx q[3];
rz(3.1252089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65288654) q[0];
sx q[0];
rz(-1.9911433) q[0];
sx q[0];
rz(2.0892129) q[0];
rz(1.9309689) q[1];
sx q[1];
rz(-1.724842) q[1];
sx q[1];
rz(0.88422424) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8540871) q[0];
sx q[0];
rz(-0.26212087) q[0];
sx q[0];
rz(-0.8192807) q[0];
rz(0.82489478) q[2];
sx q[2];
rz(-2.5755304) q[2];
sx q[2];
rz(-0.60264665) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9235232) q[1];
sx q[1];
rz(-1.4265718) q[1];
sx q[1];
rz(-1.6378535) q[1];
x q[2];
rz(2.1518425) q[3];
sx q[3];
rz(-2.5269066) q[3];
sx q[3];
rz(2.1637672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9050682) q[2];
sx q[2];
rz(-0.69513598) q[2];
sx q[2];
rz(2.0489676) q[2];
rz(0.24400273) q[3];
sx q[3];
rz(-0.87385333) q[3];
sx q[3];
rz(-2.3905335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7364175) q[0];
sx q[0];
rz(-3.1389132) q[0];
sx q[0];
rz(1.8238235) q[0];
rz(-0.70619839) q[1];
sx q[1];
rz(-1.462992) q[1];
sx q[1];
rz(-0.96907369) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9861525) q[0];
sx q[0];
rz(-1.6797721) q[0];
sx q[0];
rz(-3.100813) q[0];
rz(-pi) q[1];
rz(0.71408748) q[2];
sx q[2];
rz(-1.6220777) q[2];
sx q[2];
rz(0.72061348) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5507257) q[1];
sx q[1];
rz(-2.496965) q[1];
sx q[1];
rz(-2.4025737) q[1];
rz(-pi) q[2];
rz(-3*pi/11) q[3];
sx q[3];
rz(-2.487605) q[3];
sx q[3];
rz(1.6884782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6270854) q[2];
sx q[2];
rz(-2.5532477) q[2];
sx q[2];
rz(2.7039841) q[2];
rz(3.0833516) q[3];
sx q[3];
rz(-2.2831459) q[3];
sx q[3];
rz(-1.6102012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3095734) q[0];
sx q[0];
rz(-2.11634) q[0];
sx q[0];
rz(1.7818041) q[0];
rz(-1.4490022) q[1];
sx q[1];
rz(-2.2429008) q[1];
sx q[1];
rz(-1.70599) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5933701) q[0];
sx q[0];
rz(-2.3832294) q[0];
sx q[0];
rz(-0.81292696) q[0];
x q[1];
rz(-0.82201634) q[2];
sx q[2];
rz(-3.1146345) q[2];
sx q[2];
rz(2.6798623) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3064733) q[1];
sx q[1];
rz(-1.8715579) q[1];
sx q[1];
rz(-1.1919828) q[1];
rz(-pi) q[2];
rz(-2.2613951) q[3];
sx q[3];
rz(-0.62303632) q[3];
sx q[3];
rz(-0.33720371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.74063611) q[2];
sx q[2];
rz(-1.8813958) q[2];
sx q[2];
rz(0.17399542) q[2];
rz(-2.1334355) q[3];
sx q[3];
rz(-0.62767902) q[3];
sx q[3];
rz(-0.15667668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.85132861) q[0];
sx q[0];
rz(-0.86306089) q[0];
sx q[0];
rz(-0.1329578) q[0];
rz(-0.48775396) q[1];
sx q[1];
rz(-2.1552591) q[1];
sx q[1];
rz(-1.3652323) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2030905) q[0];
sx q[0];
rz(-1.7068958) q[0];
sx q[0];
rz(3.0654728) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6985967) q[2];
sx q[2];
rz(-1.5489849) q[2];
sx q[2];
rz(0.70805659) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.21309585) q[1];
sx q[1];
rz(-2.1017385) q[1];
sx q[1];
rz(-1.5099768) q[1];
rz(2.7216464) q[3];
sx q[3];
rz(-1.7405602) q[3];
sx q[3];
rz(1.1743634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0662971) q[2];
sx q[2];
rz(-1.2434881) q[2];
sx q[2];
rz(0.42262849) q[2];
rz(-2.1832809) q[3];
sx q[3];
rz(-3.002353) q[3];
sx q[3];
rz(0.2376093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1388824) q[0];
sx q[0];
rz(-0.30160987) q[0];
sx q[0];
rz(2.8310006) q[0];
rz(2.8839135) q[1];
sx q[1];
rz(-1.6380761) q[1];
sx q[1];
rz(-2.2176567) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63147488) q[0];
sx q[0];
rz(-0.81561136) q[0];
sx q[0];
rz(-2.5328013) q[0];
rz(-pi) q[1];
rz(1.5827492) q[2];
sx q[2];
rz(-2.0593615) q[2];
sx q[2];
rz(-1.919463) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.82211557) q[1];
sx q[1];
rz(-1.9648106) q[1];
sx q[1];
rz(1.0072717) q[1];
rz(0.71420788) q[3];
sx q[3];
rz(-2.0651157) q[3];
sx q[3];
rz(-2.2480272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6442287) q[2];
sx q[2];
rz(-2.6022537) q[2];
sx q[2];
rz(1.7473934) q[2];
rz(-1.9723069) q[3];
sx q[3];
rz(-1.0401657) q[3];
sx q[3];
rz(2.6031301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7681463) q[0];
sx q[0];
rz(-3.0420619) q[0];
sx q[0];
rz(1.3884397) q[0];
rz(-3.1346079) q[1];
sx q[1];
rz(-0.81022898) q[1];
sx q[1];
rz(3.1034234) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5725806) q[0];
sx q[0];
rz(-2.5323212) q[0];
sx q[0];
rz(-1.6174181) q[0];
x q[1];
rz(-1.1881234) q[2];
sx q[2];
rz(-1.4854991) q[2];
sx q[2];
rz(1.485699) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0109509) q[1];
sx q[1];
rz(-2.7656733) q[1];
sx q[1];
rz(-1.4569267) q[1];
rz(2.5358701) q[3];
sx q[3];
rz(-1.0014135) q[3];
sx q[3];
rz(3.0063418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.96597153) q[2];
sx q[2];
rz(-1.9922545) q[2];
sx q[2];
rz(1.1981298) q[2];
rz(-1.0839869) q[3];
sx q[3];
rz(-0.67892781) q[3];
sx q[3];
rz(-0.23588022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1643628) q[0];
sx q[0];
rz(-1.9259763) q[0];
sx q[0];
rz(1.5019793) q[0];
rz(1.4981131) q[1];
sx q[1];
rz(-0.5935946) q[1];
sx q[1];
rz(-0.53131663) q[1];
rz(-0.48508587) q[2];
sx q[2];
rz(-0.1242287) q[2];
sx q[2];
rz(2.6835174) q[2];
rz(-2.2723972) q[3];
sx q[3];
rz(-1.8164608) q[3];
sx q[3];
rz(1.7432004) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
