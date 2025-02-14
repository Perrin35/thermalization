OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1563675) q[0];
sx q[0];
rz(-1.2824143) q[0];
sx q[0];
rz(-0.13089827) q[0];
rz(0.52783293) q[1];
sx q[1];
rz(3.5117709) q[1];
sx q[1];
rz(12.389046) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41577121) q[0];
sx q[0];
rz(-2.3292589) q[0];
sx q[0];
rz(-1.9363296) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.22819569) q[2];
sx q[2];
rz(-1.8642117) q[2];
sx q[2];
rz(0.55390893) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6244391) q[1];
sx q[1];
rz(-2.6047047) q[1];
sx q[1];
rz(-1.0978903) q[1];
rz(-pi) q[2];
rz(0.99448326) q[3];
sx q[3];
rz(-1.3399933) q[3];
sx q[3];
rz(-0.9524012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6191972) q[2];
sx q[2];
rz(-1.088524) q[2];
sx q[2];
rz(1.8001451) q[2];
rz(1.3522735) q[3];
sx q[3];
rz(-1.6100581) q[3];
sx q[3];
rz(1.8488098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3905447) q[0];
sx q[0];
rz(-1.9123257) q[0];
sx q[0];
rz(-0.51502645) q[0];
rz(0.36188778) q[1];
sx q[1];
rz(-1.7444976) q[1];
sx q[1];
rz(-2.1451758) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2856355) q[0];
sx q[0];
rz(-2.8295316) q[0];
sx q[0];
rz(2.2770004) q[0];
rz(0.17536945) q[2];
sx q[2];
rz(-1.1423472) q[2];
sx q[2];
rz(-2.1074949) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7512618) q[1];
sx q[1];
rz(-1.0511025) q[1];
sx q[1];
rz(-0.10607509) q[1];
rz(0.89869638) q[3];
sx q[3];
rz(-0.65124159) q[3];
sx q[3];
rz(-0.30168331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7289875) q[2];
sx q[2];
rz(-0.94634405) q[2];
sx q[2];
rz(-1.4808572) q[2];
rz(-0.49731538) q[3];
sx q[3];
rz(-2.1931084) q[3];
sx q[3];
rz(2.4174387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5363252) q[0];
sx q[0];
rz(-1.241642) q[0];
sx q[0];
rz(-1.9507677) q[0];
rz(2.7438927) q[1];
sx q[1];
rz(-1.5813446) q[1];
sx q[1];
rz(-0.89888987) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99756344) q[0];
sx q[0];
rz(-1.9099417) q[0];
sx q[0];
rz(-3.0710286) q[0];
rz(-1.4853046) q[2];
sx q[2];
rz(-1.7608661) q[2];
sx q[2];
rz(-2.9051733) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6582322) q[1];
sx q[1];
rz(-1.0387324) q[1];
sx q[1];
rz(-0.05183713) q[1];
x q[2];
rz(1.0520027) q[3];
sx q[3];
rz(-1.4201323) q[3];
sx q[3];
rz(1.0778395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0209823) q[2];
sx q[2];
rz(-1.9526498) q[2];
sx q[2];
rz(-0.19213842) q[2];
rz(0.7297248) q[3];
sx q[3];
rz(-2.9967522) q[3];
sx q[3];
rz(-2.5016968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14908734) q[0];
sx q[0];
rz(-2.3893116) q[0];
sx q[0];
rz(1.3080904) q[0];
rz(0.84450841) q[1];
sx q[1];
rz(-1.4627855) q[1];
sx q[1];
rz(-2.5194397) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4554696) q[0];
sx q[0];
rz(-1.0884443) q[0];
sx q[0];
rz(1.1022105) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6255369) q[2];
sx q[2];
rz(-1.3856882) q[2];
sx q[2];
rz(0.4601882) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.20609186) q[1];
sx q[1];
rz(-2.969944) q[1];
sx q[1];
rz(-0.30613456) q[1];
x q[2];
rz(0.34028168) q[3];
sx q[3];
rz(-2.8064686) q[3];
sx q[3];
rz(-1.7652896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.017612351) q[2];
sx q[2];
rz(-1.7698741) q[2];
sx q[2];
rz(-2.2464216) q[2];
rz(-2.7923942) q[3];
sx q[3];
rz(-0.57217351) q[3];
sx q[3];
rz(0.23381843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(1.2772086) q[0];
sx q[0];
rz(-1.9338436) q[0];
sx q[0];
rz(-2.5982017) q[0];
rz(-0.1132938) q[1];
sx q[1];
rz(-1.3985059) q[1];
sx q[1];
rz(1.2921804) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77063507) q[0];
sx q[0];
rz(-1.1877302) q[0];
sx q[0];
rz(-2.2983453) q[0];
rz(-1.0903484) q[2];
sx q[2];
rz(-0.70864973) q[2];
sx q[2];
rz(2.6330269) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1230916) q[1];
sx q[1];
rz(-2.1827509) q[1];
sx q[1];
rz(0.83142821) q[1];
rz(-pi) q[2];
rz(-0.85815737) q[3];
sx q[3];
rz(-2.174509) q[3];
sx q[3];
rz(-1.9626372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8285404) q[2];
sx q[2];
rz(-1.7834168) q[2];
sx q[2];
rz(-0.48946112) q[2];
rz(-2.475907) q[3];
sx q[3];
rz(-2.5299215) q[3];
sx q[3];
rz(-0.78589511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2419484) q[0];
sx q[0];
rz(-0.90270942) q[0];
sx q[0];
rz(-0.06074252) q[0];
rz(1.2784866) q[1];
sx q[1];
rz(-2.2272031) q[1];
sx q[1];
rz(-2.1155105) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8908949) q[0];
sx q[0];
rz(-2.1259667) q[0];
sx q[0];
rz(-2.5722136) q[0];
rz(-pi) q[1];
rz(-0.32569484) q[2];
sx q[2];
rz(-1.0267339) q[2];
sx q[2];
rz(1.0808672) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.287331) q[1];
sx q[1];
rz(-1.320144) q[1];
sx q[1];
rz(-1.7118368) q[1];
rz(-pi) q[2];
rz(-2.6137976) q[3];
sx q[3];
rz(-1.8547117) q[3];
sx q[3];
rz(0.75263587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8817899) q[2];
sx q[2];
rz(-1.5851574) q[2];
sx q[2];
rz(-3.0628824) q[2];
rz(1.6591266) q[3];
sx q[3];
rz(-0.31646287) q[3];
sx q[3];
rz(2.6993774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92597961) q[0];
sx q[0];
rz(-2.7613566) q[0];
sx q[0];
rz(-1.6967787) q[0];
rz(2.3346057) q[1];
sx q[1];
rz(-1.3953367) q[1];
sx q[1];
rz(-0.011946202) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1769971) q[0];
sx q[0];
rz(-0.93206333) q[0];
sx q[0];
rz(-3.1082736) q[0];
x q[1];
rz(-2.6339954) q[2];
sx q[2];
rz(-1.6656808) q[2];
sx q[2];
rz(2.0213057) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.50004369) q[1];
sx q[1];
rz(-0.9608568) q[1];
sx q[1];
rz(0.088492101) q[1];
rz(1.9254382) q[3];
sx q[3];
rz(-0.4222479) q[3];
sx q[3];
rz(0.99340445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.547946) q[2];
sx q[2];
rz(-2.839489) q[2];
sx q[2];
rz(-0.83615237) q[2];
rz(-3.0892843) q[3];
sx q[3];
rz(-1.6736504) q[3];
sx q[3];
rz(-0.93594319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9713822) q[0];
sx q[0];
rz(-0.8571856) q[0];
sx q[0];
rz(-0.49825391) q[0];
rz(1.0055297) q[1];
sx q[1];
rz(-1.6906831) q[1];
sx q[1];
rz(3.0301869) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83269925) q[0];
sx q[0];
rz(-1.667262) q[0];
sx q[0];
rz(1.1497496) q[0];
x q[1];
rz(1.3206215) q[2];
sx q[2];
rz(-1.9873575) q[2];
sx q[2];
rz(-2.5377459) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1131683) q[1];
sx q[1];
rz(-2.3481124) q[1];
sx q[1];
rz(-1.9774578) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8008055) q[3];
sx q[3];
rz(-2.7772745) q[3];
sx q[3];
rz(-0.43393578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9446543) q[2];
sx q[2];
rz(-1.3613746) q[2];
sx q[2];
rz(-1.7353752) q[2];
rz(-1.1574636) q[3];
sx q[3];
rz(-2.1486053) q[3];
sx q[3];
rz(-1.5520613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.120753) q[0];
sx q[0];
rz(-0.73796219) q[0];
sx q[0];
rz(-1.4388168) q[0];
rz(-2.2976047) q[1];
sx q[1];
rz(-1.8344717) q[1];
sx q[1];
rz(-2.9356975) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5563894) q[0];
sx q[0];
rz(-2.2015338) q[0];
sx q[0];
rz(-1.6210733) q[0];
x q[1];
rz(3.0768422) q[2];
sx q[2];
rz(-2.4528385) q[2];
sx q[2];
rz(1.233686) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2515427) q[1];
sx q[1];
rz(-1.492842) q[1];
sx q[1];
rz(3.0061017) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.467942) q[3];
sx q[3];
rz(-0.60551548) q[3];
sx q[3];
rz(-0.57693627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0042808) q[2];
sx q[2];
rz(-1.0575123) q[2];
sx q[2];
rz(1.5409957) q[2];
rz(-0.48504034) q[3];
sx q[3];
rz(-1.3554327) q[3];
sx q[3];
rz(-0.11390991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.4892905) q[0];
sx q[0];
rz(-3.089383) q[0];
sx q[0];
rz(-1.8452277) q[0];
rz(-1.0507874) q[1];
sx q[1];
rz(-1.4605582) q[1];
sx q[1];
rz(-2.5349862) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25897988) q[0];
sx q[0];
rz(-1.734223) q[0];
sx q[0];
rz(3.1085186) q[0];
rz(-1.6982618) q[2];
sx q[2];
rz(-2.1268788) q[2];
sx q[2];
rz(-3.0024204) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9258825) q[1];
sx q[1];
rz(-1.5097188) q[1];
sx q[1];
rz(-0.30558773) q[1];
x q[2];
rz(-0.28429417) q[3];
sx q[3];
rz(-1.6830256) q[3];
sx q[3];
rz(-2.1865213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7217094) q[2];
sx q[2];
rz(-1.1315283) q[2];
sx q[2];
rz(-2.4427872) q[2];
rz(-1.4612259) q[3];
sx q[3];
rz(-2.1278087) q[3];
sx q[3];
rz(-0.47718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.2496495) q[0];
sx q[0];
rz(-1.5105381) q[0];
sx q[0];
rz(1.426209) q[0];
rz(-2.7632948) q[1];
sx q[1];
rz(-2.5143647) q[1];
sx q[1];
rz(-0.73077269) q[1];
rz(2.8332491) q[2];
sx q[2];
rz(-2.0258122) q[2];
sx q[2];
rz(-3.0292055) q[2];
rz(-1.5398239) q[3];
sx q[3];
rz(-2.1270558) q[3];
sx q[3];
rz(1.1598641) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
