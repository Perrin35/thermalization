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
rz(3.0313015) q[0];
sx q[0];
rz(-1.8561441) q[0];
sx q[0];
rz(0.31874803) q[0];
rz(-1.9593852) q[1];
sx q[1];
rz(-1.477834) q[1];
sx q[1];
rz(1.1629265) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0025151) q[0];
sx q[0];
rz(-2.7292508) q[0];
sx q[0];
rz(0.95100944) q[0];
rz(1.5854344) q[2];
sx q[2];
rz(-1.7649289) q[2];
sx q[2];
rz(-0.090198719) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9388814) q[1];
sx q[1];
rz(-1.2699281) q[1];
sx q[1];
rz(-1.2510514) q[1];
rz(-2.5666994) q[3];
sx q[3];
rz(-1.5539317) q[3];
sx q[3];
rz(0.3080388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.17243871) q[2];
sx q[2];
rz(-0.46010751) q[2];
sx q[2];
rz(-3.0058506) q[2];
rz(2.8067449) q[3];
sx q[3];
rz(-1.9183153) q[3];
sx q[3];
rz(2.7443583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6927476) q[0];
sx q[0];
rz(-1.0053585) q[0];
sx q[0];
rz(0.94648615) q[0];
rz(1.954156) q[1];
sx q[1];
rz(-0.58381909) q[1];
sx q[1];
rz(2.8576287) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89881247) q[0];
sx q[0];
rz(-1.3889342) q[0];
sx q[0];
rz(-0.19354958) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87209629) q[2];
sx q[2];
rz(-1.3882033) q[2];
sx q[2];
rz(-0.77699772) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.52580035) q[1];
sx q[1];
rz(-0.4471356) q[1];
sx q[1];
rz(2.8195803) q[1];
x q[2];
rz(-0.3881298) q[3];
sx q[3];
rz(-1.8947269) q[3];
sx q[3];
rz(0.38979724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9339319) q[2];
sx q[2];
rz(-1.8457103) q[2];
sx q[2];
rz(0.80061039) q[2];
rz(-1.8396395) q[3];
sx q[3];
rz(-1.1336361) q[3];
sx q[3];
rz(3.083526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0170853) q[0];
sx q[0];
rz(-0.43211102) q[0];
sx q[0];
rz(-2.3486163) q[0];
rz(2.4555581) q[1];
sx q[1];
rz(-1.425309) q[1];
sx q[1];
rz(-0.7652024) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8898123) q[0];
sx q[0];
rz(-1.3597825) q[0];
sx q[0];
rz(3.0821256) q[0];
rz(0.051570895) q[2];
sx q[2];
rz(-1.602939) q[2];
sx q[2];
rz(2.3598537) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.4884035) q[1];
sx q[1];
rz(-0.81522664) q[1];
sx q[1];
rz(1.1865739) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92242494) q[3];
sx q[3];
rz(-0.88705813) q[3];
sx q[3];
rz(2.2278663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.56696314) q[2];
sx q[2];
rz(-0.69090635) q[2];
sx q[2];
rz(-1.8951269) q[2];
rz(-1.9555107) q[3];
sx q[3];
rz(-1.7639152) q[3];
sx q[3];
rz(-0.20868364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3935811) q[0];
sx q[0];
rz(-1.3205386) q[0];
sx q[0];
rz(3.0082974) q[0];
rz(0.26946274) q[1];
sx q[1];
rz(-2.2270484) q[1];
sx q[1];
rz(-2.6752313) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9892557) q[0];
sx q[0];
rz(-1.6902802) q[0];
sx q[0];
rz(0.59520041) q[0];
x q[1];
rz(1.3106636) q[2];
sx q[2];
rz(-0.64412737) q[2];
sx q[2];
rz(0.82714236) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.129648) q[1];
sx q[1];
rz(-0.62800558) q[1];
sx q[1];
rz(2.8824214) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5138054) q[3];
sx q[3];
rz(-1.5132959) q[3];
sx q[3];
rz(-2.1138885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93126297) q[2];
sx q[2];
rz(-0.42079058) q[2];
sx q[2];
rz(-0.24410625) q[2];
rz(1.7804451) q[3];
sx q[3];
rz(-2.1972392) q[3];
sx q[3];
rz(1.2066427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72531438) q[0];
sx q[0];
rz(-2.2940574) q[0];
sx q[0];
rz(3.0552979) q[0];
rz(1.7591954) q[1];
sx q[1];
rz(-1.0603797) q[1];
sx q[1];
rz(-1.8550526) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4327398) q[0];
sx q[0];
rz(-1.4012432) q[0];
sx q[0];
rz(-1.1772593) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7930512) q[2];
sx q[2];
rz(-0.80446595) q[2];
sx q[2];
rz(-0.90057238) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2225688) q[1];
sx q[1];
rz(-1.6529539) q[1];
sx q[1];
rz(-0.8687353) q[1];
x q[2];
rz(0.87950893) q[3];
sx q[3];
rz(-1.9909166) q[3];
sx q[3];
rz(-2.0685982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7354108) q[2];
sx q[2];
rz(-1.7291131) q[2];
sx q[2];
rz(-1.1253051) q[2];
rz(-0.79648894) q[3];
sx q[3];
rz(-2.4060566) q[3];
sx q[3];
rz(-2.9430732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35561246) q[0];
sx q[0];
rz(-2.8481843) q[0];
sx q[0];
rz(-0.32817131) q[0];
rz(2.7525821) q[1];
sx q[1];
rz(-1.5704472) q[1];
sx q[1];
rz(1.6860115) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7195014) q[0];
sx q[0];
rz(-1.6515284) q[0];
sx q[0];
rz(0.21105612) q[0];
rz(-pi) q[1];
rz(3.0148618) q[2];
sx q[2];
rz(-1.7046844) q[2];
sx q[2];
rz(0.19370475) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.751669) q[1];
sx q[1];
rz(-1.8018541) q[1];
sx q[1];
rz(-0.91737813) q[1];
x q[2];
rz(-1.6809627) q[3];
sx q[3];
rz(-0.36073437) q[3];
sx q[3];
rz(-2.9025789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.97633156) q[2];
sx q[2];
rz(-2.2165522) q[2];
sx q[2];
rz(0.46249214) q[2];
rz(2.1180604) q[3];
sx q[3];
rz(-2.0868802) q[3];
sx q[3];
rz(2.3033477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8031215) q[0];
sx q[0];
rz(-1.4305038) q[0];
sx q[0];
rz(-2.9679003) q[0];
rz(2.9202785) q[1];
sx q[1];
rz(-2.4559655) q[1];
sx q[1];
rz(-1.7783222) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0039788469) q[0];
sx q[0];
rz(-0.44671808) q[0];
sx q[0];
rz(1.1398619) q[0];
rz(0.82011804) q[2];
sx q[2];
rz(-3.0228399) q[2];
sx q[2];
rz(1.135716) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1044019) q[1];
sx q[1];
rz(-2.1774182) q[1];
sx q[1];
rz(2.1527075) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2990427) q[3];
sx q[3];
rz(-1.3370965) q[3];
sx q[3];
rz(-0.2971479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.26019105) q[2];
sx q[2];
rz(-1.1611725) q[2];
sx q[2];
rz(-1.4455522) q[2];
rz(-0.77406231) q[3];
sx q[3];
rz(-2.4249707) q[3];
sx q[3];
rz(-1.6707481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.562302) q[0];
sx q[0];
rz(-0.93920541) q[0];
sx q[0];
rz(-2.8939409) q[0];
rz(-0.66894764) q[1];
sx q[1];
rz(-1.1890143) q[1];
sx q[1];
rz(0.98562366) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7108817) q[0];
sx q[0];
rz(-2.1257002) q[0];
sx q[0];
rz(1.4355833) q[0];
rz(-pi) q[1];
rz(0.47565461) q[2];
sx q[2];
rz(-0.61810571) q[2];
sx q[2];
rz(2.3007768) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0541773) q[1];
sx q[1];
rz(-1.2654516) q[1];
sx q[1];
rz(-1.1609689) q[1];
rz(-pi) q[2];
rz(-0.59496112) q[3];
sx q[3];
rz(-1.5719126) q[3];
sx q[3];
rz(3.0240455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0011657) q[2];
sx q[2];
rz(-0.75158921) q[2];
sx q[2];
rz(-1.2933732) q[2];
rz(1.2398531) q[3];
sx q[3];
rz(-0.69609061) q[3];
sx q[3];
rz(-2.4333439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2769315) q[0];
sx q[0];
rz(-1.7129352) q[0];
sx q[0];
rz(2.1760333) q[0];
rz(0.37636617) q[1];
sx q[1];
rz(-1.8195567) q[1];
sx q[1];
rz(-1.9115062) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9133643) q[0];
sx q[0];
rz(-1.3051864) q[0];
sx q[0];
rz(1.5188541) q[0];
rz(-1.240991) q[2];
sx q[2];
rz(-1.9422741) q[2];
sx q[2];
rz(-0.35201752) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4415746) q[1];
sx q[1];
rz(-2.5315365) q[1];
sx q[1];
rz(3.0185591) q[1];
rz(-3.0885124) q[3];
sx q[3];
rz(-2.3276668) q[3];
sx q[3];
rz(2.2102714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2824715) q[2];
sx q[2];
rz(-0.66994795) q[2];
sx q[2];
rz(0.13151375) q[2];
rz(2.6604743) q[3];
sx q[3];
rz(-2.2115464) q[3];
sx q[3];
rz(-0.49829811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-1.1199101) q[0];
sx q[0];
rz(-0.57761884) q[0];
sx q[0];
rz(0.48900327) q[0];
rz(1.8148212) q[1];
sx q[1];
rz(-1.2811456) q[1];
sx q[1];
rz(2.9723523) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45825935) q[0];
sx q[0];
rz(-1.1348327) q[0];
sx q[0];
rz(-3.0100214) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6851114) q[2];
sx q[2];
rz(-2.940753) q[2];
sx q[2];
rz(1.072509) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3551363) q[1];
sx q[1];
rz(-0.82477714) q[1];
sx q[1];
rz(2.5591947) q[1];
rz(-pi) q[2];
rz(-2.333913) q[3];
sx q[3];
rz(-0.39356222) q[3];
sx q[3];
rz(0.53787947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5532316) q[2];
sx q[2];
rz(-2.0560122) q[2];
sx q[2];
rz(-2.9300743) q[2];
rz(1.5254947) q[3];
sx q[3];
rz(-2.1414521) q[3];
sx q[3];
rz(-0.26743993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78032988) q[0];
sx q[0];
rz(-1.4365256) q[0];
sx q[0];
rz(-2.8509675) q[0];
rz(0.81193874) q[1];
sx q[1];
rz(-0.62807905) q[1];
sx q[1];
rz(-1.5608578) q[1];
rz(1.590325) q[2];
sx q[2];
rz(-1.1371798) q[2];
sx q[2];
rz(-0.22424998) q[2];
rz(-0.62150443) q[3];
sx q[3];
rz(-0.92172289) q[3];
sx q[3];
rz(-0.48529101) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
