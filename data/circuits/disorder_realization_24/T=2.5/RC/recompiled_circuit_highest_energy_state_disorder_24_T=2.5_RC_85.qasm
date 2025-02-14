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
rz(2.6296122) q[0];
sx q[0];
rz(-0.57096243) q[0];
sx q[0];
rz(-0.1877187) q[0];
rz(-2.3277148) q[1];
sx q[1];
rz(-1.4717646) q[1];
sx q[1];
rz(-1.2893113) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16351249) q[0];
sx q[0];
rz(-2.8239282) q[0];
sx q[0];
rz(-2.4207522) q[0];
x q[1];
rz(-3.0953636) q[2];
sx q[2];
rz(-0.57114375) q[2];
sx q[2];
rz(2.2276304) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.047886176) q[1];
sx q[1];
rz(-1.1125065) q[1];
sx q[1];
rz(1.8480193) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1407152) q[3];
sx q[3];
rz(-2.7288247) q[3];
sx q[3];
rz(1.4980396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2671555) q[2];
sx q[2];
rz(-2.0384553) q[2];
sx q[2];
rz(0.77679408) q[2];
rz(-1.8393501) q[3];
sx q[3];
rz(-1.4812508) q[3];
sx q[3];
rz(1.2031901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60748196) q[0];
sx q[0];
rz(-2.6597839) q[0];
sx q[0];
rz(2.1433461) q[0];
rz(2.7718995) q[1];
sx q[1];
rz(-1.0414711) q[1];
sx q[1];
rz(-2.0236156) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6551483) q[0];
sx q[0];
rz(-1.0308497) q[0];
sx q[0];
rz(0.47294954) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68685617) q[2];
sx q[2];
rz(-1.1026898) q[2];
sx q[2];
rz(1.8091701) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0740023) q[1];
sx q[1];
rz(-2.6704881) q[1];
sx q[1];
rz(-1.8995238) q[1];
rz(-0.26882986) q[3];
sx q[3];
rz(-2.1850366) q[3];
sx q[3];
rz(1.5419095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.86576858) q[2];
sx q[2];
rz(-1.7663225) q[2];
sx q[2];
rz(-2.7562874) q[2];
rz(-3.1028808) q[3];
sx q[3];
rz(-2.7888515) q[3];
sx q[3];
rz(-2.501131) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63967079) q[0];
sx q[0];
rz(-2.6646035) q[0];
sx q[0];
rz(-1.84024) q[0];
rz(3.0838857) q[1];
sx q[1];
rz(-2.6951908) q[1];
sx q[1];
rz(1.3267964) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0327137) q[0];
sx q[0];
rz(-1.098806) q[0];
sx q[0];
rz(-0.95979877) q[0];
rz(2.4270664) q[2];
sx q[2];
rz(-0.21103141) q[2];
sx q[2];
rz(2.2843366) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4918265) q[1];
sx q[1];
rz(-2.5297102) q[1];
sx q[1];
rz(-0.51458451) q[1];
rz(-pi) q[2];
rz(-2.2928724) q[3];
sx q[3];
rz(-1.8406788) q[3];
sx q[3];
rz(0.11972846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3441299) q[2];
sx q[2];
rz(-0.8351438) q[2];
sx q[2];
rz(-0.75378913) q[2];
rz(1.7720743) q[3];
sx q[3];
rz(-1.6794208) q[3];
sx q[3];
rz(0.65521017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9264483) q[0];
sx q[0];
rz(-1.0855874) q[0];
sx q[0];
rz(1.7468859) q[0];
rz(1.9649547) q[1];
sx q[1];
rz(-1.6329012) q[1];
sx q[1];
rz(-0.36697695) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1934051) q[0];
sx q[0];
rz(-0.087554878) q[0];
sx q[0];
rz(2.2746536) q[0];
rz(-2.8072678) q[2];
sx q[2];
rz(-1.8211289) q[2];
sx q[2];
rz(-1.3368397) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.61287294) q[1];
sx q[1];
rz(-2.7320128) q[1];
sx q[1];
rz(-2.6367841) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6567201) q[3];
sx q[3];
rz(-1.20487) q[3];
sx q[3];
rz(-2.576862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7299812) q[2];
sx q[2];
rz(-1.9133762) q[2];
sx q[2];
rz(1.8709987) q[2];
rz(-0.96108428) q[3];
sx q[3];
rz(-1.7226487) q[3];
sx q[3];
rz(0.050475033) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5807895) q[0];
sx q[0];
rz(-2.2152948) q[0];
sx q[0];
rz(1.3128989) q[0];
rz(1.1543697) q[1];
sx q[1];
rz(-1.1416953) q[1];
sx q[1];
rz(-2.5837574) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2656888) q[0];
sx q[0];
rz(-0.77123986) q[0];
sx q[0];
rz(-0.79637595) q[0];
rz(-0.16061546) q[2];
sx q[2];
rz(-1.6470419) q[2];
sx q[2];
rz(0.74197021) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.085943) q[1];
sx q[1];
rz(-0.80971566) q[1];
sx q[1];
rz(-1.147406) q[1];
x q[2];
rz(-2.6799503) q[3];
sx q[3];
rz(-1.7064693) q[3];
sx q[3];
rz(-1.5380579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4981726) q[2];
sx q[2];
rz(-1.8489445) q[2];
sx q[2];
rz(2.2922929) q[2];
rz(0.15527209) q[3];
sx q[3];
rz(-0.97584358) q[3];
sx q[3];
rz(0.37080216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.3575386) q[0];
sx q[0];
rz(-0.58670601) q[0];
sx q[0];
rz(0.46911711) q[0];
rz(0.47438374) q[1];
sx q[1];
rz(-0.7862888) q[1];
sx q[1];
rz(0.75538409) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4159214) q[0];
sx q[0];
rz(-1.5752107) q[0];
sx q[0];
rz(1.2995385) q[0];
x q[1];
rz(0.62894507) q[2];
sx q[2];
rz(-1.8417995) q[2];
sx q[2];
rz(2.6586201) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50690097) q[1];
sx q[1];
rz(-1.697487) q[1];
sx q[1];
rz(-1.3562528) q[1];
rz(-2.7128025) q[3];
sx q[3];
rz(-2.2615823) q[3];
sx q[3];
rz(-0.59861983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0834121) q[2];
sx q[2];
rz(-1.8184793) q[2];
sx q[2];
rz(0.18784909) q[2];
rz(-1.8958873) q[3];
sx q[3];
rz(-1.7626423) q[3];
sx q[3];
rz(1.0904306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(1.93194) q[0];
sx q[0];
rz(-1.8745475) q[0];
sx q[0];
rz(0.39091045) q[0];
rz(0.93005013) q[1];
sx q[1];
rz(-1.4314194) q[1];
sx q[1];
rz(-1.8375058) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.025607312) q[0];
sx q[0];
rz(-1.5918333) q[0];
sx q[0];
rz(-1.9217092) q[0];
x q[1];
rz(-1.485059) q[2];
sx q[2];
rz(-1.9112327) q[2];
sx q[2];
rz(0.96735937) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3871349) q[1];
sx q[1];
rz(-1.4031193) q[1];
sx q[1];
rz(-1.1134336) q[1];
rz(-pi) q[2];
rz(3.0797187) q[3];
sx q[3];
rz(-0.29682595) q[3];
sx q[3];
rz(-2.5421531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8730674) q[2];
sx q[2];
rz(-2.8748942) q[2];
sx q[2];
rz(-3.0282057) q[2];
rz(0.015965613) q[3];
sx q[3];
rz(-1.6458052) q[3];
sx q[3];
rz(-0.16460831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3575344) q[0];
sx q[0];
rz(-1.4736195) q[0];
sx q[0];
rz(-3.0175324) q[0];
rz(-0.1217753) q[1];
sx q[1];
rz(-2.3901794) q[1];
sx q[1];
rz(1.6990936) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7989726) q[0];
sx q[0];
rz(-0.39726394) q[0];
sx q[0];
rz(1.8333866) q[0];
rz(-pi) q[1];
rz(0.14367468) q[2];
sx q[2];
rz(-0.29479237) q[2];
sx q[2];
rz(-0.93570342) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.21181078) q[1];
sx q[1];
rz(-0.36961296) q[1];
sx q[1];
rz(-1.9464284) q[1];
rz(0.27713953) q[3];
sx q[3];
rz(-2.3921514) q[3];
sx q[3];
rz(-2.998874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1759935) q[2];
sx q[2];
rz(-1.3600574) q[2];
sx q[2];
rz(0.74481258) q[2];
rz(-2.2143769) q[3];
sx q[3];
rz(-1.5085446) q[3];
sx q[3];
rz(2.0492699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6532779) q[0];
sx q[0];
rz(-1.7019685) q[0];
sx q[0];
rz(-0.22536817) q[0];
rz(2.0291746) q[1];
sx q[1];
rz(-0.77443361) q[1];
sx q[1];
rz(-0.89881277) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3254125) q[0];
sx q[0];
rz(-1.0409779) q[0];
sx q[0];
rz(-2.2384032) q[0];
rz(-pi) q[1];
rz(-0.65943879) q[2];
sx q[2];
rz(-2.5213679) q[2];
sx q[2];
rz(2.8218215) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0045709) q[1];
sx q[1];
rz(-2.4413476) q[1];
sx q[1];
rz(-1.2546468) q[1];
rz(-pi) q[2];
rz(1.5360654) q[3];
sx q[3];
rz(-2.8485549) q[3];
sx q[3];
rz(1.422783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0972458) q[2];
sx q[2];
rz(-1.4298507) q[2];
sx q[2];
rz(-0.35935768) q[2];
rz(2.8906631) q[3];
sx q[3];
rz(-1.072262) q[3];
sx q[3];
rz(0.76771626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3971685) q[0];
sx q[0];
rz(-3.011062) q[0];
sx q[0];
rz(-2.2644444) q[0];
rz(1.5646704) q[1];
sx q[1];
rz(-1.6404057) q[1];
sx q[1];
rz(2.3559779) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025185275) q[0];
sx q[0];
rz(-1.8078363) q[0];
sx q[0];
rz(-2.9171506) q[0];
rz(-2.6076786) q[2];
sx q[2];
rz(-1.1964239) q[2];
sx q[2];
rz(-1.716734) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4686376) q[1];
sx q[1];
rz(-1.5002999) q[1];
sx q[1];
rz(-0.34480178) q[1];
x q[2];
rz(-0.50650017) q[3];
sx q[3];
rz(-1.439038) q[3];
sx q[3];
rz(-1.7817287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2239573) q[2];
sx q[2];
rz(-0.79006299) q[2];
sx q[2];
rz(-2.4033974) q[2];
rz(-2.2186642) q[3];
sx q[3];
rz(-1.8332558) q[3];
sx q[3];
rz(-0.2963399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.7947163) q[0];
sx q[0];
rz(-1.7807757) q[0];
sx q[0];
rz(-2.1697252) q[0];
rz(-0.90732668) q[1];
sx q[1];
rz(-2.0569888) q[1];
sx q[1];
rz(2.8203698) q[1];
rz(-0.39225929) q[2];
sx q[2];
rz(-1.0259368) q[2];
sx q[2];
rz(-2.1403014) q[2];
rz(2.3774556) q[3];
sx q[3];
rz(-2.4485179) q[3];
sx q[3];
rz(-2.5701523) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
