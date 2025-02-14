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
rz(1.6990868) q[0];
sx q[0];
rz(-1.8449755) q[0];
sx q[0];
rz(-1.698864) q[0];
rz(-0.35248414) q[1];
sx q[1];
rz(-0.52630693) q[1];
sx q[1];
rz(-1.3679282) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099771899) q[0];
sx q[0];
rz(-1.213014) q[0];
sx q[0];
rz(1.8463184) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6085195) q[2];
sx q[2];
rz(-1.2153271) q[2];
sx q[2];
rz(1.0743574) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.753073) q[1];
sx q[1];
rz(-2.6021829) q[1];
sx q[1];
rz(2.8864157) q[1];
x q[2];
rz(-0.97729965) q[3];
sx q[3];
rz(-2.233089) q[3];
sx q[3];
rz(-2.4031201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0736531) q[2];
sx q[2];
rz(-2.0906134) q[2];
sx q[2];
rz(-1.2037207) q[2];
rz(-2.1667571) q[3];
sx q[3];
rz(-0.9674415) q[3];
sx q[3];
rz(1.2459285) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3881653) q[0];
sx q[0];
rz(-2.0797256) q[0];
sx q[0];
rz(-0.88428307) q[0];
rz(2.4321804) q[1];
sx q[1];
rz(-1.246289) q[1];
sx q[1];
rz(2.2394004) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7912238) q[0];
sx q[0];
rz(-2.02105) q[0];
sx q[0];
rz(1.3421935) q[0];
rz(3.0386417) q[2];
sx q[2];
rz(-0.96625096) q[2];
sx q[2];
rz(1.1510731) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1556405) q[1];
sx q[1];
rz(-1.2725755) q[1];
sx q[1];
rz(-1.31485) q[1];
rz(1.7092) q[3];
sx q[3];
rz(-1.4908199) q[3];
sx q[3];
rz(-3.0113335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0238637) q[2];
sx q[2];
rz(-0.29568299) q[2];
sx q[2];
rz(-1.4884865) q[2];
rz(1.9272517) q[3];
sx q[3];
rz(-1.6818654) q[3];
sx q[3];
rz(1.2795718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0073256) q[0];
sx q[0];
rz(-1.7420344) q[0];
sx q[0];
rz(2.9752327) q[0];
rz(0.69794377) q[1];
sx q[1];
rz(-0.78015399) q[1];
sx q[1];
rz(1.4770329) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0255942) q[0];
sx q[0];
rz(-2.6995772) q[0];
sx q[0];
rz(1.5777977) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8516305) q[2];
sx q[2];
rz(-1.3094433) q[2];
sx q[2];
rz(-0.42119831) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.229847) q[1];
sx q[1];
rz(-2.1087136) q[1];
sx q[1];
rz(2.6078546) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6432041) q[3];
sx q[3];
rz(-1.5880463) q[3];
sx q[3];
rz(-1.1422865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.6537689) q[2];
sx q[2];
rz(-1.1710125) q[2];
sx q[2];
rz(-2.5269395) q[2];
rz(1.4608308) q[3];
sx q[3];
rz(-1.5771644) q[3];
sx q[3];
rz(-3.0998668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6093269) q[0];
sx q[0];
rz(-1.8121239) q[0];
sx q[0];
rz(-1.6229269) q[0];
rz(-0.80865771) q[1];
sx q[1];
rz(-1.2963908) q[1];
sx q[1];
rz(-3.0928968) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8762465) q[0];
sx q[0];
rz(-2.2314203) q[0];
sx q[0];
rz(1.1424078) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28566177) q[2];
sx q[2];
rz(-2.7451519) q[2];
sx q[2];
rz(-0.91914058) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.21203707) q[1];
sx q[1];
rz(-1.0921282) q[1];
sx q[1];
rz(-0.57600601) q[1];
x q[2];
rz(-1.0297431) q[3];
sx q[3];
rz(-0.67922878) q[3];
sx q[3];
rz(-1.7558869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.101717) q[2];
sx q[2];
rz(-0.81284916) q[2];
sx q[2];
rz(1.8854878) q[2];
rz(3.0806372) q[3];
sx q[3];
rz(-2.239949) q[3];
sx q[3];
rz(-0.92756334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1561279) q[0];
sx q[0];
rz(-1.1277072) q[0];
sx q[0];
rz(0.10966478) q[0];
rz(-2.1127286) q[1];
sx q[1];
rz(-1.2867462) q[1];
sx q[1];
rz(-1.5922155) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4745195) q[0];
sx q[0];
rz(-1.1231849) q[0];
sx q[0];
rz(1.2569129) q[0];
rz(-pi) q[1];
rz(-2.778901) q[2];
sx q[2];
rz(-1.0981993) q[2];
sx q[2];
rz(-2.9031799) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.60051892) q[1];
sx q[1];
rz(-1.1364577) q[1];
sx q[1];
rz(0.51323607) q[1];
x q[2];
rz(-3.1072633) q[3];
sx q[3];
rz(-1.5331556) q[3];
sx q[3];
rz(-3.056789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5900383) q[2];
sx q[2];
rz(-0.2912713) q[2];
sx q[2];
rz(0.48198286) q[2];
rz(-0.41989741) q[3];
sx q[3];
rz(-1.0651383) q[3];
sx q[3];
rz(-0.70986748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9919306) q[0];
sx q[0];
rz(-0.88352942) q[0];
sx q[0];
rz(-2.4290207) q[0];
rz(0.86992162) q[1];
sx q[1];
rz(-2.7674119) q[1];
sx q[1];
rz(3.1178927) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7256692) q[0];
sx q[0];
rz(-2.5953889) q[0];
sx q[0];
rz(0.31809455) q[0];
rz(-pi) q[1];
rz(-0.65448241) q[2];
sx q[2];
rz(-0.81771446) q[2];
sx q[2];
rz(0.47374642) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.99147532) q[1];
sx q[1];
rz(-1.5722701) q[1];
sx q[1];
rz(1.7365843) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38857821) q[3];
sx q[3];
rz(-2.6286308) q[3];
sx q[3];
rz(1.4298317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3370543) q[2];
sx q[2];
rz(-2.4961175) q[2];
sx q[2];
rz(-2.0013334) q[2];
rz(-1.1537457) q[3];
sx q[3];
rz(-1.9269582) q[3];
sx q[3];
rz(-1.8203189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-1.7158647) q[0];
sx q[0];
rz(-1.3874929) q[0];
sx q[0];
rz(-2.7737889) q[0];
rz(-1.4999207) q[1];
sx q[1];
rz(-2.2190614) q[1];
sx q[1];
rz(1.0949162) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.80847) q[0];
sx q[0];
rz(-0.54352353) q[0];
sx q[0];
rz(1.1103866) q[0];
rz(-pi) q[1];
rz(2.5911683) q[2];
sx q[2];
rz(-1.4604521) q[2];
sx q[2];
rz(0.36729022) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3125632) q[1];
sx q[1];
rz(-2.9051648) q[1];
sx q[1];
rz(-0.076024012) q[1];
x q[2];
rz(-1.6206018) q[3];
sx q[3];
rz(-0.40400716) q[3];
sx q[3];
rz(-1.2329519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0423923) q[2];
sx q[2];
rz(-1.3513869) q[2];
sx q[2];
rz(-3.0858827) q[2];
rz(0.67692155) q[3];
sx q[3];
rz(-0.61546314) q[3];
sx q[3];
rz(-2.2630579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.95213503) q[0];
sx q[0];
rz(-0.5624693) q[0];
sx q[0];
rz(-2.1761555) q[0];
rz(1.2646487) q[1];
sx q[1];
rz(-2.3955884) q[1];
sx q[1];
rz(-0.3269349) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1345987) q[0];
sx q[0];
rz(-2.7732458) q[0];
sx q[0];
rz(-1.8828431) q[0];
rz(2.6956396) q[2];
sx q[2];
rz(-1.5805681) q[2];
sx q[2];
rz(-0.31644145) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2911243) q[1];
sx q[1];
rz(-1.7479436) q[1];
sx q[1];
rz(-1.3590823) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.502874) q[3];
sx q[3];
rz(-2.0188031) q[3];
sx q[3];
rz(0.94326708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38291976) q[2];
sx q[2];
rz(-0.70432538) q[2];
sx q[2];
rz(0.79603377) q[2];
rz(-0.087513611) q[3];
sx q[3];
rz(-0.60860601) q[3];
sx q[3];
rz(2.2039738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4336808) q[0];
sx q[0];
rz(-1.3734564) q[0];
sx q[0];
rz(0.32980907) q[0];
rz(3.0682796) q[1];
sx q[1];
rz(-2.4344756) q[1];
sx q[1];
rz(2.0475533) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46991321) q[0];
sx q[0];
rz(-1.6080583) q[0];
sx q[0];
rz(0.47084634) q[0];
rz(-pi) q[1];
rz(-1.173063) q[2];
sx q[2];
rz(-1.5088044) q[2];
sx q[2];
rz(-2.6336489) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0713404) q[1];
sx q[1];
rz(-1.0254745) q[1];
sx q[1];
rz(-1.4163744) q[1];
x q[2];
rz(1.1393896) q[3];
sx q[3];
rz(-1.6069222) q[3];
sx q[3];
rz(-2.2808035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0426992) q[2];
sx q[2];
rz(-2.3697479) q[2];
sx q[2];
rz(-1.137255) q[2];
rz(0.62816652) q[3];
sx q[3];
rz(-0.38574949) q[3];
sx q[3];
rz(-1.0073193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6657418) q[0];
sx q[0];
rz(-2.6306212) q[0];
sx q[0];
rz(1.092859) q[0];
rz(-1.8765556) q[1];
sx q[1];
rz(-1.3053514) q[1];
sx q[1];
rz(-2.1078033) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2726934) q[0];
sx q[0];
rz(-1.0995691) q[0];
sx q[0];
rz(-0.97288218) q[0];
x q[1];
rz(-0.25569465) q[2];
sx q[2];
rz(-0.77633023) q[2];
sx q[2];
rz(1.2344373) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.30822291) q[1];
sx q[1];
rz(-2.3607681) q[1];
sx q[1];
rz(0.60365795) q[1];
rz(1.3218811) q[3];
sx q[3];
rz(-2.9985709) q[3];
sx q[3];
rz(1.420457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8826302) q[2];
sx q[2];
rz(-0.7236824) q[2];
sx q[2];
rz(2.9222729) q[2];
rz(-0.80260459) q[3];
sx q[3];
rz(-1.31253) q[3];
sx q[3];
rz(-2.820211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5103067) q[0];
sx q[0];
rz(-1.318537) q[0];
sx q[0];
rz(1.0986811) q[0];
rz(2.9639099) q[1];
sx q[1];
rz(-1.0164574) q[1];
sx q[1];
rz(-0.16998092) q[1];
rz(3.0776824) q[2];
sx q[2];
rz(-0.8693246) q[2];
sx q[2];
rz(-1.0988591) q[2];
rz(-2.1228409) q[3];
sx q[3];
rz(-0.22291796) q[3];
sx q[3];
rz(2.1414584) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
