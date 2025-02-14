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
rz(0.91020838) q[0];
sx q[0];
rz(-0.75054032) q[0];
sx q[0];
rz(-2.5810177) q[0];
rz(1.9380467) q[1];
sx q[1];
rz(0.37625852) q[1];
sx q[1];
rz(9.6179554) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7734186) q[0];
sx q[0];
rz(-0.76245284) q[0];
sx q[0];
rz(-0.66526757) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7142263) q[2];
sx q[2];
rz(-0.71038112) q[2];
sx q[2];
rz(1.8153035) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5467599) q[1];
sx q[1];
rz(-1.7103737) q[1];
sx q[1];
rz(-0.91800331) q[1];
rz(0.47115342) q[3];
sx q[3];
rz(-0.36590016) q[3];
sx q[3];
rz(0.19972502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.91118497) q[2];
sx q[2];
rz(-2.513803) q[2];
sx q[2];
rz(-0.5578624) q[2];
rz(1.2281536) q[3];
sx q[3];
rz(-1.6817254) q[3];
sx q[3];
rz(0.094999878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-1.2896344) q[0];
sx q[0];
rz(-1.8486706) q[0];
sx q[0];
rz(1.5648382) q[0];
rz(-1.1467038) q[1];
sx q[1];
rz(-1.7161938) q[1];
sx q[1];
rz(-2.1099405) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1250336) q[0];
sx q[0];
rz(-1.3368589) q[0];
sx q[0];
rz(1.0256605) q[0];
rz(2.7286058) q[2];
sx q[2];
rz(-1.7888165) q[2];
sx q[2];
rz(-0.08139164) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37687846) q[1];
sx q[1];
rz(-0.8974613) q[1];
sx q[1];
rz(-2.7002579) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50526039) q[3];
sx q[3];
rz(-1.7499515) q[3];
sx q[3];
rz(0.0091088692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8138294) q[2];
sx q[2];
rz(-0.29080614) q[2];
sx q[2];
rz(-1.7247058) q[2];
rz(-2.318577) q[3];
sx q[3];
rz(-1.6133285) q[3];
sx q[3];
rz(-2.4970162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3642035) q[0];
sx q[0];
rz(-2.0222029) q[0];
sx q[0];
rz(-0.35225824) q[0];
rz(0.38905713) q[1];
sx q[1];
rz(-1.9720826) q[1];
sx q[1];
rz(-0.22638098) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8733559) q[0];
sx q[0];
rz(-1.653329) q[0];
sx q[0];
rz(-1.8186452) q[0];
rz(-pi) q[1];
rz(2.0820322) q[2];
sx q[2];
rz(-0.63768688) q[2];
sx q[2];
rz(-2.9342655) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.375294) q[1];
sx q[1];
rz(-2.0336478) q[1];
sx q[1];
rz(-2.0867324) q[1];
rz(1.4949293) q[3];
sx q[3];
rz(-2.0056723) q[3];
sx q[3];
rz(-2.7135389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9809197) q[2];
sx q[2];
rz(-2.1143165) q[2];
sx q[2];
rz(2.7799907) q[2];
rz(2.8412039) q[3];
sx q[3];
rz(-1.8301679) q[3];
sx q[3];
rz(-0.96295199) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2196197) q[0];
sx q[0];
rz(-1.0924871) q[0];
sx q[0];
rz(0.12119448) q[0];
rz(-1.9425862) q[1];
sx q[1];
rz(-1.0795178) q[1];
sx q[1];
rz(-1.2746864) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29783937) q[0];
sx q[0];
rz(-0.22704253) q[0];
sx q[0];
rz(0.20926306) q[0];
rz(-pi) q[1];
rz(-0.077353296) q[2];
sx q[2];
rz(-2.935545) q[2];
sx q[2];
rz(-1.8872583) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9642942) q[1];
sx q[1];
rz(-1.7907447) q[1];
sx q[1];
rz(0.58348685) q[1];
rz(-1.2043456) q[3];
sx q[3];
rz(-2.7525558) q[3];
sx q[3];
rz(-1.6477444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.049383) q[2];
sx q[2];
rz(-1.4468687) q[2];
sx q[2];
rz(1.3775187) q[2];
rz(-2.0130017) q[3];
sx q[3];
rz(-2.1994574) q[3];
sx q[3];
rz(0.82730627) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89151299) q[0];
sx q[0];
rz(-0.87123195) q[0];
sx q[0];
rz(-1.0846035) q[0];
rz(0.67101038) q[1];
sx q[1];
rz(-1.6295461) q[1];
sx q[1];
rz(-0.074140851) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0947572) q[0];
sx q[0];
rz(-2.1291385) q[0];
sx q[0];
rz(1.895491) q[0];
x q[1];
rz(-2.8507892) q[2];
sx q[2];
rz(-1.8755696) q[2];
sx q[2];
rz(1.0195315) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3931294) q[1];
sx q[1];
rz(-0.65915758) q[1];
sx q[1];
rz(1.1200957) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1187333) q[3];
sx q[3];
rz(-0.30395711) q[3];
sx q[3];
rz(0.077465103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.025617754) q[2];
sx q[2];
rz(-0.74883777) q[2];
sx q[2];
rz(-2.1785114) q[2];
rz(-1.7668308) q[3];
sx q[3];
rz(-1.129496) q[3];
sx q[3];
rz(0.53205427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9328203) q[0];
sx q[0];
rz(-2.8498579) q[0];
sx q[0];
rz(1.2579086) q[0];
rz(0.84016291) q[1];
sx q[1];
rz(-1.9636619) q[1];
sx q[1];
rz(-2.449583) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35273472) q[0];
sx q[0];
rz(-1.2458572) q[0];
sx q[0];
rz(0.28577013) q[0];
x q[1];
rz(1.7501159) q[2];
sx q[2];
rz(-0.81294927) q[2];
sx q[2];
rz(-1.6008962) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.51541182) q[1];
sx q[1];
rz(-1.0728991) q[1];
sx q[1];
rz(0.84705686) q[1];
rz(0.70051381) q[3];
sx q[3];
rz(-0.76863241) q[3];
sx q[3];
rz(-0.045846102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.96364) q[2];
sx q[2];
rz(-1.0963564) q[2];
sx q[2];
rz(-0.5160416) q[2];
rz(-1.3456723) q[3];
sx q[3];
rz(-2.7005152) q[3];
sx q[3];
rz(-2.1399982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26821414) q[0];
sx q[0];
rz(-0.28609797) q[0];
sx q[0];
rz(1.0299261) q[0];
rz(-2.4870807) q[1];
sx q[1];
rz(-1.3348568) q[1];
sx q[1];
rz(-2.9676504) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3315863) q[0];
sx q[0];
rz(-0.86649202) q[0];
sx q[0];
rz(0.84789168) q[0];
rz(-pi) q[1];
rz(1.0368749) q[2];
sx q[2];
rz(-0.838111) q[2];
sx q[2];
rz(1.378841) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.084981266) q[1];
sx q[1];
rz(-2.1826996) q[1];
sx q[1];
rz(2.8065744) q[1];
rz(-1.8344363) q[3];
sx q[3];
rz(-1.0705433) q[3];
sx q[3];
rz(-1.5280346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8339463) q[2];
sx q[2];
rz(-2.6070194) q[2];
sx q[2];
rz(-1.3391116) q[2];
rz(0.80859679) q[3];
sx q[3];
rz(-1.1491038) q[3];
sx q[3];
rz(1.8554525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22063743) q[0];
sx q[0];
rz(-1.0458825) q[0];
sx q[0];
rz(-1.9761696) q[0];
rz(2.9479345) q[1];
sx q[1];
rz(-0.8332738) q[1];
sx q[1];
rz(-1.9906893) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1068368) q[0];
sx q[0];
rz(-2.2118705) q[0];
sx q[0];
rz(1.8888372) q[0];
x q[1];
rz(1.3940349) q[2];
sx q[2];
rz(-2.3675248) q[2];
sx q[2];
rz(1.7947527) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.47025679) q[1];
sx q[1];
rz(-2.0681948) q[1];
sx q[1];
rz(-0.48723826) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64784785) q[3];
sx q[3];
rz(-1.0937249) q[3];
sx q[3];
rz(1.939524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.2151486) q[2];
sx q[2];
rz(-1.0215267) q[2];
sx q[2];
rz(-1.8683757) q[2];
rz(2.6269954) q[3];
sx q[3];
rz(-0.319258) q[3];
sx q[3];
rz(0.74015051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1505245) q[0];
sx q[0];
rz(-0.06572289) q[0];
sx q[0];
rz(0.42732987) q[0];
rz(-1.7006251) q[1];
sx q[1];
rz(-1.7040375) q[1];
sx q[1];
rz(0.89868054) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5367765) q[0];
sx q[0];
rz(-1.9990078) q[0];
sx q[0];
rz(-2.8003119) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9814453) q[2];
sx q[2];
rz(-1.169983) q[2];
sx q[2];
rz(-3.1360195) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8391621) q[1];
sx q[1];
rz(-1.3962455) q[1];
sx q[1];
rz(-1.3329092) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4472404) q[3];
sx q[3];
rz(-2.4789841) q[3];
sx q[3];
rz(2.724805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23540363) q[2];
sx q[2];
rz(-2.0960505) q[2];
sx q[2];
rz(1.40353) q[2];
rz(-0.67052001) q[3];
sx q[3];
rz(-1.5454005) q[3];
sx q[3];
rz(-1.8286573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8197166) q[0];
sx q[0];
rz(-2.1724367) q[0];
sx q[0];
rz(-2.4153391) q[0];
rz(2.4655474) q[1];
sx q[1];
rz(-1.36422) q[1];
sx q[1];
rz(-0.77686754) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1438734) q[0];
sx q[0];
rz(-1.6164427) q[0];
sx q[0];
rz(-2.2111228) q[0];
rz(-1.5258342) q[2];
sx q[2];
rz(-1.4248214) q[2];
sx q[2];
rz(-3.1342497) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90112359) q[1];
sx q[1];
rz(-1.8778342) q[1];
sx q[1];
rz(-1.9528051) q[1];
rz(-pi) q[2];
rz(-2.5489574) q[3];
sx q[3];
rz(-2.0505954) q[3];
sx q[3];
rz(2.0247839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1399416) q[2];
sx q[2];
rz(-2.9652014) q[2];
sx q[2];
rz(-1.6528992) q[2];
rz(-1.1267003) q[3];
sx q[3];
rz(-1.1688346) q[3];
sx q[3];
rz(-2.6039629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5269932) q[0];
sx q[0];
rz(-0.64170964) q[0];
sx q[0];
rz(1.5668305) q[0];
rz(-2.3552786) q[1];
sx q[1];
rz(-2.96824) q[1];
sx q[1];
rz(2.7460964) q[1];
rz(0.98914374) q[2];
sx q[2];
rz(-2.157515) q[2];
sx q[2];
rz(-0.6872057) q[2];
rz(2.8948404) q[3];
sx q[3];
rz(-2.3507531) q[3];
sx q[3];
rz(1.7274461) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
