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
rz(0.93208575) q[0];
sx q[0];
rz(4.1756364) q[0];
sx q[0];
rz(11.812165) q[0];
rz(1.7475313) q[1];
sx q[1];
rz(-0.59352195) q[1];
sx q[1];
rz(-1.4374179) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93255723) q[0];
sx q[0];
rz(-0.13929312) q[0];
sx q[0];
rz(-0.4362696) q[0];
rz(-1.6637592) q[2];
sx q[2];
rz(-0.39080253) q[2];
sx q[2];
rz(1.2338232) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61753786) q[1];
sx q[1];
rz(-2.1270942) q[1];
sx q[1];
rz(1.2134882) q[1];
rz(-pi) q[2];
rz(1.2487605) q[3];
sx q[3];
rz(-0.40169558) q[3];
sx q[3];
rz(-0.63689771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6393911) q[2];
sx q[2];
rz(-2.907967) q[2];
sx q[2];
rz(2.8761253) q[2];
rz(1.4207077) q[3];
sx q[3];
rz(-1.9749494) q[3];
sx q[3];
rz(1.407912) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0625286) q[0];
sx q[0];
rz(-1.205235) q[0];
sx q[0];
rz(-2.089654) q[0];
rz(2.1417292) q[1];
sx q[1];
rz(-1.5971767) q[1];
sx q[1];
rz(2.3725407) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5239983) q[0];
sx q[0];
rz(-2.4193561) q[0];
sx q[0];
rz(-2.2666032) q[0];
x q[1];
rz(0.090754358) q[2];
sx q[2];
rz(-0.57951515) q[2];
sx q[2];
rz(-1.7837634) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5271142) q[1];
sx q[1];
rz(-1.1243781) q[1];
sx q[1];
rz(-2.9698257) q[1];
rz(-pi) q[2];
rz(-2.6399986) q[3];
sx q[3];
rz(-2.6644649) q[3];
sx q[3];
rz(0.17233822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.96180463) q[2];
sx q[2];
rz(-2.3369117) q[2];
sx q[2];
rz(-1.556832) q[2];
rz(0.70476091) q[3];
sx q[3];
rz(-2.9289991) q[3];
sx q[3];
rz(1.0889277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8162808) q[0];
sx q[0];
rz(-0.75958696) q[0];
sx q[0];
rz(-2.3749206) q[0];
rz(2.7491772) q[1];
sx q[1];
rz(-0.86154834) q[1];
sx q[1];
rz(-0.79889417) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93358675) q[0];
sx q[0];
rz(-2.0547935) q[0];
sx q[0];
rz(0.010864929) q[0];
rz(0.25071745) q[2];
sx q[2];
rz(-1.483416) q[2];
sx q[2];
rz(-2.6911497) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0943999) q[1];
sx q[1];
rz(-2.5421738) q[1];
sx q[1];
rz(2.1381738) q[1];
rz(-pi) q[2];
rz(1.464099) q[3];
sx q[3];
rz(-1.7472526) q[3];
sx q[3];
rz(-0.7430232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.40156349) q[2];
sx q[2];
rz(-2.8758958) q[2];
sx q[2];
rz(3.1166039) q[2];
rz(-1.7449024) q[3];
sx q[3];
rz(-1.2324421) q[3];
sx q[3];
rz(-0.11740824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5609189) q[0];
sx q[0];
rz(-2.6011401) q[0];
sx q[0];
rz(-0.29676357) q[0];
rz(-0.98776039) q[1];
sx q[1];
rz(-1.8901653) q[1];
sx q[1];
rz(-0.74772778) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0714617) q[0];
sx q[0];
rz(-0.80118766) q[0];
sx q[0];
rz(1.2642994) q[0];
rz(0.50570391) q[2];
sx q[2];
rz(-0.78185588) q[2];
sx q[2];
rz(3.0145558) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8750413) q[1];
sx q[1];
rz(-1.5593441) q[1];
sx q[1];
rz(-2.7391866) q[1];
rz(-pi) q[2];
rz(-2.1652664) q[3];
sx q[3];
rz(-2.1552999) q[3];
sx q[3];
rz(2.2849831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1838386) q[2];
sx q[2];
rz(-0.87205333) q[2];
sx q[2];
rz(-2.6217065) q[2];
rz(0.94684354) q[3];
sx q[3];
rz(-1.1917453) q[3];
sx q[3];
rz(0.051518353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2648322) q[0];
sx q[0];
rz(-2.2417289) q[0];
sx q[0];
rz(1.1967891) q[0];
rz(-1.4747249) q[1];
sx q[1];
rz(-2.3366172) q[1];
sx q[1];
rz(2.8428452) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0550514) q[0];
sx q[0];
rz(-2.3712681) q[0];
sx q[0];
rz(2.6280759) q[0];
x q[1];
rz(0.49017723) q[2];
sx q[2];
rz(-0.85974795) q[2];
sx q[2];
rz(-2.4012964) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1890002) q[1];
sx q[1];
rz(-2.3680284) q[1];
sx q[1];
rz(-0.1632593) q[1];
rz(-2.0386054) q[3];
sx q[3];
rz(-0.75689471) q[3];
sx q[3];
rz(-1.5278969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3195191) q[2];
sx q[2];
rz(-0.75054979) q[2];
sx q[2];
rz(0.015241148) q[2];
rz(-3.0139253) q[3];
sx q[3];
rz(-0.36104194) q[3];
sx q[3];
rz(-2.291919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78758883) q[0];
sx q[0];
rz(-0.51743999) q[0];
sx q[0];
rz(-2.2678243) q[0];
rz(-1.9173701) q[1];
sx q[1];
rz(-1.0226378) q[1];
sx q[1];
rz(1.9220985) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3412832) q[0];
sx q[0];
rz(-0.62200802) q[0];
sx q[0];
rz(2.3013297) q[0];
rz(-2.4151426) q[2];
sx q[2];
rz(-1.1442351) q[2];
sx q[2];
rz(-2.5202765) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1186824) q[1];
sx q[1];
rz(-1.7798276) q[1];
sx q[1];
rz(0.19281705) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55762724) q[3];
sx q[3];
rz(-2.4708797) q[3];
sx q[3];
rz(2.022555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.445861) q[2];
sx q[2];
rz(-1.6709238) q[2];
sx q[2];
rz(-0.50103465) q[2];
rz(-1.3387574) q[3];
sx q[3];
rz(-2.7300291) q[3];
sx q[3];
rz(-1.9921654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5434791) q[0];
sx q[0];
rz(-1.9444332) q[0];
sx q[0];
rz(0.57394779) q[0];
rz(2.5698938) q[1];
sx q[1];
rz(-2.1015002) q[1];
sx q[1];
rz(-1.5843102) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1645419) q[0];
sx q[0];
rz(-0.91616154) q[0];
sx q[0];
rz(-0.47619168) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1102067) q[2];
sx q[2];
rz(-0.71515036) q[2];
sx q[2];
rz(-0.88727027) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7148278) q[1];
sx q[1];
rz(-1.0958156) q[1];
sx q[1];
rz(-2.2248041) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9628223) q[3];
sx q[3];
rz(-1.9147263) q[3];
sx q[3];
rz(1.0052538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.456936) q[2];
sx q[2];
rz(-1.9921649) q[2];
sx q[2];
rz(1.1535025) q[2];
rz(-1.5270816) q[3];
sx q[3];
rz(-1.0228415) q[3];
sx q[3];
rz(1.3326741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1077147) q[0];
sx q[0];
rz(-0.77924538) q[0];
sx q[0];
rz(-0.23767924) q[0];
rz(2.4029845) q[1];
sx q[1];
rz(-1.9269201) q[1];
sx q[1];
rz(3.1225263) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6929853) q[0];
sx q[0];
rz(-0.073090471) q[0];
sx q[0];
rz(-1.684171) q[0];
rz(-pi) q[1];
rz(-2.3735061) q[2];
sx q[2];
rz(-1.0407018) q[2];
sx q[2];
rz(-1.4370611) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.42603675) q[1];
sx q[1];
rz(-1.5320679) q[1];
sx q[1];
rz(-1.4891796) q[1];
x q[2];
rz(-0.98293368) q[3];
sx q[3];
rz(-1.0140099) q[3];
sx q[3];
rz(2.9012321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1261403) q[2];
sx q[2];
rz(-1.7917874) q[2];
sx q[2];
rz(-0.17507412) q[2];
rz(-1.3779047) q[3];
sx q[3];
rz(-2.521069) q[3];
sx q[3];
rz(-0.11848816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0835251) q[0];
sx q[0];
rz(-1.7836934) q[0];
sx q[0];
rz(0.51496664) q[0];
rz(0.99634755) q[1];
sx q[1];
rz(-2.0774272) q[1];
sx q[1];
rz(-1.6641585) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9436059) q[0];
sx q[0];
rz(-1.0062402) q[0];
sx q[0];
rz(-2.415231) q[0];
rz(-pi) q[1];
rz(0.049749497) q[2];
sx q[2];
rz(-2.5526415) q[2];
sx q[2];
rz(-2.4204554) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1783289) q[1];
sx q[1];
rz(-2.43578) q[1];
sx q[1];
rz(-2.3030445) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9733125) q[3];
sx q[3];
rz(-2.5096143) q[3];
sx q[3];
rz(-3.0807854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0763863) q[2];
sx q[2];
rz(-1.4896769) q[2];
sx q[2];
rz(-0.15929407) q[2];
rz(1.6431036) q[3];
sx q[3];
rz(-2.4703333) q[3];
sx q[3];
rz(1.338965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0945187) q[0];
sx q[0];
rz(-3.0975332) q[0];
sx q[0];
rz(0.86167589) q[0];
rz(1.2791951) q[1];
sx q[1];
rz(-2.8586614) q[1];
sx q[1];
rz(-0.18712015) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1094799) q[0];
sx q[0];
rz(-0.36046991) q[0];
sx q[0];
rz(-0.69860639) q[0];
rz(-pi) q[1];
rz(1.2995903) q[2];
sx q[2];
rz(-1.091963) q[2];
sx q[2];
rz(1.9128654) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5622065) q[1];
sx q[1];
rz(-0.77695266) q[1];
sx q[1];
rz(-0.38511606) q[1];
rz(-1.2544823) q[3];
sx q[3];
rz(-1.1939002) q[3];
sx q[3];
rz(1.0939456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6992135) q[2];
sx q[2];
rz(-2.1670161) q[2];
sx q[2];
rz(1.7787735) q[2];
rz(1.6393225) q[3];
sx q[3];
rz(-2.5876741) q[3];
sx q[3];
rz(-1.7207918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5300071) q[0];
sx q[0];
rz(-1.8066318) q[0];
sx q[0];
rz(-3.1351177) q[0];
rz(0.50637983) q[1];
sx q[1];
rz(-0.19229278) q[1];
sx q[1];
rz(-2.2813588) q[1];
rz(1.8851595) q[2];
sx q[2];
rz(-2.6994575) q[2];
sx q[2];
rz(-0.14002249) q[2];
rz(-2.8741638) q[3];
sx q[3];
rz(-1.3602644) q[3];
sx q[3];
rz(-0.2841831) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
