OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.78046507) q[0];
sx q[0];
rz(3.830885) q[0];
sx q[0];
rz(9.7552714) q[0];
rz(0.32980546) q[1];
sx q[1];
rz(-0.84996119) q[1];
sx q[1];
rz(0.70911521) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.37492) q[0];
sx q[0];
rz(-1.2776889) q[0];
sx q[0];
rz(0.93107443) q[0];
rz(-pi) q[1];
rz(1.66967) q[2];
sx q[2];
rz(-2.8266202) q[2];
sx q[2];
rz(1.0613943) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4830556) q[1];
sx q[1];
rz(-0.72209789) q[1];
sx q[1];
rz(0.87941054) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37681864) q[3];
sx q[3];
rz(-2.0599277) q[3];
sx q[3];
rz(0.34987846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4101397) q[2];
sx q[2];
rz(-0.4099161) q[2];
sx q[2];
rz(1.5343792) q[2];
rz(2.2065227) q[3];
sx q[3];
rz(-1.8362703) q[3];
sx q[3];
rz(0.7888166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7222897) q[0];
sx q[0];
rz(-0.062148217) q[0];
sx q[0];
rz(0.62227917) q[0];
rz(0.17624804) q[1];
sx q[1];
rz(-1.9269678) q[1];
sx q[1];
rz(0.91631779) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8371671) q[0];
sx q[0];
rz(-2.2080748) q[0];
sx q[0];
rz(-0.62563719) q[0];
rz(-pi) q[1];
rz(1.4257405) q[2];
sx q[2];
rz(-0.74160355) q[2];
sx q[2];
rz(1.2765826) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2434477) q[1];
sx q[1];
rz(-1.8971649) q[1];
sx q[1];
rz(1.9279724) q[1];
rz(1.6244435) q[3];
sx q[3];
rz(-2.7430153) q[3];
sx q[3];
rz(0.22565354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.24094412) q[2];
sx q[2];
rz(-0.99664656) q[2];
sx q[2];
rz(-0.82143482) q[2];
rz(3.1243096) q[3];
sx q[3];
rz(-0.79764962) q[3];
sx q[3];
rz(0.78330529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.9538552) q[0];
sx q[0];
rz(-1.5783577) q[0];
sx q[0];
rz(2.1333372) q[0];
rz(3.1058274) q[1];
sx q[1];
rz(-1.4859896) q[1];
sx q[1];
rz(0.52454138) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13428155) q[0];
sx q[0];
rz(-1.6005922) q[0];
sx q[0];
rz(-1.7344463) q[0];
x q[1];
rz(-1.8345941) q[2];
sx q[2];
rz(-0.909415) q[2];
sx q[2];
rz(-1.7922572) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.2136038) q[1];
sx q[1];
rz(-0.065918006) q[1];
sx q[1];
rz(-0.35857486) q[1];
rz(-pi) q[2];
rz(1.9906304) q[3];
sx q[3];
rz(-1.435278) q[3];
sx q[3];
rz(-2.6672222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3699469) q[2];
sx q[2];
rz(-1.5180072) q[2];
sx q[2];
rz(1.4952205) q[2];
rz(1.3211936) q[3];
sx q[3];
rz(-1.4097872) q[3];
sx q[3];
rz(-0.43740073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8111073) q[0];
sx q[0];
rz(-0.91902584) q[0];
sx q[0];
rz(0.088949732) q[0];
rz(-0.51070172) q[1];
sx q[1];
rz(-0.82534868) q[1];
sx q[1];
rz(-2.451992) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3084429) q[0];
sx q[0];
rz(-1.2763378) q[0];
sx q[0];
rz(-1.740728) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7961411) q[2];
sx q[2];
rz(-1.0190522) q[2];
sx q[2];
rz(-1.4622886) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.888962) q[1];
sx q[1];
rz(-1.8969715) q[1];
sx q[1];
rz(-0.77146448) q[1];
rz(-pi) q[2];
rz(-1.7161937) q[3];
sx q[3];
rz(-1.0905915) q[3];
sx q[3];
rz(-0.2916382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.66578635) q[2];
sx q[2];
rz(-0.67725724) q[2];
sx q[2];
rz(-0.62292567) q[2];
rz(2.0056491) q[3];
sx q[3];
rz(-0.1853075) q[3];
sx q[3];
rz(-2.6749271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2899807) q[0];
sx q[0];
rz(-1.8988134) q[0];
sx q[0];
rz(-0.583453) q[0];
rz(1.9955697) q[1];
sx q[1];
rz(-1.4910411) q[1];
sx q[1];
rz(-1.6437644) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0935055) q[0];
sx q[0];
rz(-1.830528) q[0];
sx q[0];
rz(0.78608677) q[0];
rz(-pi) q[1];
rz(0.096502467) q[2];
sx q[2];
rz(-0.8000024) q[2];
sx q[2];
rz(2.2068791) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6240546) q[1];
sx q[1];
rz(-2.2741286) q[1];
sx q[1];
rz(1.6066949) q[1];
rz(-pi) q[2];
rz(-1.4086401) q[3];
sx q[3];
rz(-2.080653) q[3];
sx q[3];
rz(-1.7061403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4862165) q[2];
sx q[2];
rz(-1.5950173) q[2];
sx q[2];
rz(-1.5636469) q[2];
rz(-0.90562138) q[3];
sx q[3];
rz(-0.27035299) q[3];
sx q[3];
rz(2.5031228) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49801302) q[0];
sx q[0];
rz(-1.6828515) q[0];
sx q[0];
rz(0.046982732) q[0];
rz(-0.14818305) q[1];
sx q[1];
rz(-0.89887416) q[1];
sx q[1];
rz(-1.4354338) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5934138) q[0];
sx q[0];
rz(-2.3313064) q[0];
sx q[0];
rz(2.9368648) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0211421) q[2];
sx q[2];
rz(-1.6209941) q[2];
sx q[2];
rz(-2.2113138) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1286436) q[1];
sx q[1];
rz(-1.9777021) q[1];
sx q[1];
rz(3.0055771) q[1];
x q[2];
rz(1.3330323) q[3];
sx q[3];
rz(-1.5974345) q[3];
sx q[3];
rz(2.3029072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0662213) q[2];
sx q[2];
rz(-2.1263289) q[2];
sx q[2];
rz(3.0701239) q[2];
rz(-1.6890769) q[3];
sx q[3];
rz(-1.6966597) q[3];
sx q[3];
rz(-2.215109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(2.8554095) q[0];
sx q[0];
rz(-1.8137285) q[0];
sx q[0];
rz(-0.57762161) q[0];
rz(-1.2795992) q[1];
sx q[1];
rz(-1.7956693) q[1];
sx q[1];
rz(2.1320027) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1560695) q[0];
sx q[0];
rz(-2.0445604) q[0];
sx q[0];
rz(1.9457293) q[0];
rz(-pi) q[1];
rz(1.3426443) q[2];
sx q[2];
rz(-1.9462799) q[2];
sx q[2];
rz(-1.0588888) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2958888) q[1];
sx q[1];
rz(-1.0755952) q[1];
sx q[1];
rz(1.1843029) q[1];
x q[2];
rz(0.76241242) q[3];
sx q[3];
rz(-2.0250642) q[3];
sx q[3];
rz(-1.7319958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0916831) q[2];
sx q[2];
rz(-0.35353264) q[2];
sx q[2];
rz(1.1996777) q[2];
rz(0.47232929) q[3];
sx q[3];
rz(-1.4940558) q[3];
sx q[3];
rz(-2.1388617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49884477) q[0];
sx q[0];
rz(-1.6413178) q[0];
sx q[0];
rz(2.902466) q[0];
rz(2.3587976) q[1];
sx q[1];
rz(-2.6233964) q[1];
sx q[1];
rz(2.696864) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77138222) q[0];
sx q[0];
rz(-1.0078197) q[0];
sx q[0];
rz(2.6315106) q[0];
x q[1];
rz(0.93034805) q[2];
sx q[2];
rz(-2.1313416) q[2];
sx q[2];
rz(-2.408574) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4152894) q[1];
sx q[1];
rz(-2.2518034) q[1];
sx q[1];
rz(-2.9773832) q[1];
rz(2.8430311) q[3];
sx q[3];
rz(-2.7136554) q[3];
sx q[3];
rz(0.096979389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.42177054) q[2];
sx q[2];
rz(-0.82227102) q[2];
sx q[2];
rz(-2.3262809) q[2];
rz(2.1067965) q[3];
sx q[3];
rz(-1.5650322) q[3];
sx q[3];
rz(-1.8243272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8274882) q[0];
sx q[0];
rz(-1.0359456) q[0];
sx q[0];
rz(-0.28717336) q[0];
rz(-2.9526967) q[1];
sx q[1];
rz(-2.4177528) q[1];
sx q[1];
rz(-0.33219355) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4591601) q[0];
sx q[0];
rz(-2.3258665) q[0];
sx q[0];
rz(-1.1343603) q[0];
x q[1];
rz(3.095093) q[2];
sx q[2];
rz(-0.87247712) q[2];
sx q[2];
rz(2.7547714) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0041973) q[1];
sx q[1];
rz(-1.157837) q[1];
sx q[1];
rz(-1.2910299) q[1];
rz(0.53211777) q[3];
sx q[3];
rz(-1.2174774) q[3];
sx q[3];
rz(1.4766828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2087848) q[2];
sx q[2];
rz(-2.1366426) q[2];
sx q[2];
rz(-2.3020111) q[2];
rz(1.2906637) q[3];
sx q[3];
rz(-1.3495812) q[3];
sx q[3];
rz(2.4850142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0861417) q[0];
sx q[0];
rz(-0.76776004) q[0];
sx q[0];
rz(2.0822051) q[0];
rz(-2.1620031) q[1];
sx q[1];
rz(-1.1971808) q[1];
sx q[1];
rz(-0.25451452) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5352288) q[0];
sx q[0];
rz(-1.4974721) q[0];
sx q[0];
rz(1.801977) q[0];
rz(-pi) q[1];
rz(2.7024686) q[2];
sx q[2];
rz(-0.52999485) q[2];
sx q[2];
rz(0.32967552) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6470681) q[1];
sx q[1];
rz(-1.8956603) q[1];
sx q[1];
rz(1.6560133) q[1];
rz(-pi) q[2];
rz(-0.20930807) q[3];
sx q[3];
rz(-2.6946687) q[3];
sx q[3];
rz(-0.79093864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0467726) q[2];
sx q[2];
rz(-2.5929055) q[2];
sx q[2];
rz(1.0478896) q[2];
rz(1.3607599) q[3];
sx q[3];
rz(-0.69793099) q[3];
sx q[3];
rz(0.60539436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.027325252) q[0];
sx q[0];
rz(-1.1189168) q[0];
sx q[0];
rz(3.0615321) q[0];
rz(-0.36021532) q[1];
sx q[1];
rz(-1.4604912) q[1];
sx q[1];
rz(2.1866658) q[1];
rz(1.1709471) q[2];
sx q[2];
rz(-2.5794537) q[2];
sx q[2];
rz(-0.07515547) q[2];
rz(-1.8923106) q[3];
sx q[3];
rz(-0.63890639) q[3];
sx q[3];
rz(-1.3756868) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
