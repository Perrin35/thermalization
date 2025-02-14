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
rz(0.56057492) q[0];
rz(-1.203546) q[1];
sx q[1];
rz(-0.37625852) q[1];
sx q[1];
rz(-2.9484152) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3681741) q[0];
sx q[0];
rz(-0.76245284) q[0];
sx q[0];
rz(2.4763251) q[0];
x q[1];
rz(3.019252) q[2];
sx q[2];
rz(-2.2723778) q[2];
sx q[2];
rz(2.0035494) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.59483278) q[1];
sx q[1];
rz(-1.7103737) q[1];
sx q[1];
rz(-0.91800331) q[1];
rz(-pi) q[2];
rz(-2.812593) q[3];
sx q[3];
rz(-1.4076678) q[3];
sx q[3];
rz(1.3265058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.91118497) q[2];
sx q[2];
rz(-2.513803) q[2];
sx q[2];
rz(2.5837303) q[2];
rz(1.2281536) q[3];
sx q[3];
rz(-1.4598673) q[3];
sx q[3];
rz(3.0465928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8519583) q[0];
sx q[0];
rz(-1.8486706) q[0];
sx q[0];
rz(-1.5767545) q[0];
rz(-1.1467038) q[1];
sx q[1];
rz(-1.4253989) q[1];
sx q[1];
rz(2.1099405) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0165591) q[0];
sx q[0];
rz(-1.8047338) q[0];
sx q[0];
rz(1.0256605) q[0];
rz(-pi) q[1];
rz(-0.41298683) q[2];
sx q[2];
rz(-1.7888165) q[2];
sx q[2];
rz(3.060201) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.116392) q[1];
sx q[1];
rz(-2.355819) q[1];
sx q[1];
rz(-2.0624195) q[1];
rz(-pi) q[2];
rz(-0.50526039) q[3];
sx q[3];
rz(-1.3916411) q[3];
sx q[3];
rz(-0.0091088692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.32776323) q[2];
sx q[2];
rz(-0.29080614) q[2];
sx q[2];
rz(-1.7247058) q[2];
rz(2.318577) q[3];
sx q[3];
rz(-1.6133285) q[3];
sx q[3];
rz(-0.64457646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3642035) q[0];
sx q[0];
rz(-2.0222029) q[0];
sx q[0];
rz(0.35225824) q[0];
rz(-2.7525355) q[1];
sx q[1];
rz(-1.16951) q[1];
sx q[1];
rz(-2.9152117) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3234179) q[0];
sx q[0];
rz(-1.8177839) q[0];
sx q[0];
rz(0.085121827) q[0];
x q[1];
rz(0.99708544) q[2];
sx q[2];
rz(-1.8663532) q[2];
sx q[2];
rz(-1.7869345) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5851262) q[1];
sx q[1];
rz(-2.0279998) q[1];
sx q[1];
rz(0.5208421) q[1];
rz(1.4949293) q[3];
sx q[3];
rz(-1.1359204) q[3];
sx q[3];
rz(-0.42805373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.160673) q[2];
sx q[2];
rz(-1.0272762) q[2];
sx q[2];
rz(2.7799907) q[2];
rz(-2.8412039) q[3];
sx q[3];
rz(-1.3114248) q[3];
sx q[3];
rz(2.1786407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2196197) q[0];
sx q[0];
rz(-2.0491056) q[0];
sx q[0];
rz(0.12119448) q[0];
rz(-1.1990064) q[1];
sx q[1];
rz(-2.0620748) q[1];
sx q[1];
rz(1.8669063) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.629144) q[0];
sx q[0];
rz(-1.7928018) q[0];
sx q[0];
rz(-1.52284) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0642394) q[2];
sx q[2];
rz(-2.935545) q[2];
sx q[2];
rz(-1.2543343) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9642942) q[1];
sx q[1];
rz(-1.3508479) q[1];
sx q[1];
rz(-0.58348685) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.995756) q[3];
sx q[3];
rz(-1.2088299) q[3];
sx q[3];
rz(-1.2545619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.092209665) q[2];
sx q[2];
rz(-1.694724) q[2];
sx q[2];
rz(-1.3775187) q[2];
rz(1.1285909) q[3];
sx q[3];
rz(-0.94213525) q[3];
sx q[3];
rz(2.3142864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2500797) q[0];
sx q[0];
rz(-2.2703607) q[0];
sx q[0];
rz(-2.0569892) q[0];
rz(0.67101038) q[1];
sx q[1];
rz(-1.6295461) q[1];
sx q[1];
rz(3.0674518) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6607802) q[0];
sx q[0];
rz(-0.63711053) q[0];
sx q[0];
rz(2.6694032) q[0];
rz(-pi) q[1];
rz(1.2535278) q[2];
sx q[2];
rz(-1.2937577) q[2];
sx q[2];
rz(-2.679897) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.843815) q[1];
sx q[1];
rz(-0.98688302) q[1];
sx q[1];
rz(0.32547863) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8458141) q[3];
sx q[3];
rz(-1.4396808) q[3];
sx q[3];
rz(1.0594378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.025617754) q[2];
sx q[2];
rz(-2.3927549) q[2];
sx q[2];
rz(2.1785114) q[2];
rz(-1.3747619) q[3];
sx q[3];
rz(-2.0120967) q[3];
sx q[3];
rz(-2.6095384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.2087723) q[0];
sx q[0];
rz(-2.8498579) q[0];
sx q[0];
rz(-1.2579086) q[0];
rz(2.3014297) q[1];
sx q[1];
rz(-1.9636619) q[1];
sx q[1];
rz(2.449583) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0449311) q[0];
sx q[0];
rz(-2.7122444) q[0];
sx q[0];
rz(2.2675687) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18628405) q[2];
sx q[2];
rz(-2.3669503) q[2];
sx q[2];
rz(-1.8587405) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.51541182) q[1];
sx q[1];
rz(-1.0728991) q[1];
sx q[1];
rz(2.2945358) q[1];
rz(1.0133842) q[3];
sx q[3];
rz(-1.0104826) q[3];
sx q[3];
rz(-0.91065613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.96364) q[2];
sx q[2];
rz(-2.0452363) q[2];
sx q[2];
rz(-2.625551) q[2];
rz(1.7959203) q[3];
sx q[3];
rz(-0.44107744) q[3];
sx q[3];
rz(-1.0015944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26821414) q[0];
sx q[0];
rz(-0.28609797) q[0];
sx q[0];
rz(1.0299261) q[0];
rz(-0.65451199) q[1];
sx q[1];
rz(-1.8067358) q[1];
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
rz(-0.51474606) q[2];
sx q[2];
rz(-0.87640773) q[2];
sx q[2];
rz(1.0389164) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.084981266) q[1];
sx q[1];
rz(-2.1826996) q[1];
sx q[1];
rz(0.33501825) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44486041) q[3];
sx q[3];
rz(-2.5813817) q[3];
sx q[3];
rz(-2.1261393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.30764636) q[2];
sx q[2];
rz(-0.53457326) q[2];
sx q[2];
rz(-1.3391116) q[2];
rz(-0.80859679) q[3];
sx q[3];
rz(-1.1491038) q[3];
sx q[3];
rz(1.2861402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22063743) q[0];
sx q[0];
rz(-2.0957102) q[0];
sx q[0];
rz(-1.165423) q[0];
rz(-0.19365817) q[1];
sx q[1];
rz(-2.3083189) q[1];
sx q[1];
rz(-1.1509034) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1068368) q[0];
sx q[0];
rz(-0.92972219) q[0];
sx q[0];
rz(-1.8888372) q[0];
rz(-pi) q[1];
rz(-1.3940349) q[2];
sx q[2];
rz(-0.77406787) q[2];
sx q[2];
rz(-1.34684) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.833646) q[1];
sx q[1];
rz(-0.6816136) q[1];
sx q[1];
rz(0.8591842) q[1];
rz(-pi) q[2];
rz(0.70827534) q[3];
sx q[3];
rz(-2.357956) q[3];
sx q[3];
rz(2.2275138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.2151486) q[2];
sx q[2];
rz(-2.1200659) q[2];
sx q[2];
rz(1.8683757) q[2];
rz(-2.6269954) q[3];
sx q[3];
rz(-0.319258) q[3];
sx q[3];
rz(-0.74015051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99106818) q[0];
sx q[0];
rz(-0.06572289) q[0];
sx q[0];
rz(0.42732987) q[0];
rz(1.4409675) q[1];
sx q[1];
rz(-1.7040375) q[1];
sx q[1];
rz(-2.2429121) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81954992) q[0];
sx q[0];
rz(-1.8801831) q[0];
sx q[0];
rz(1.1197107) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1601474) q[2];
sx q[2];
rz(-1.9716096) q[2];
sx q[2];
rz(-0.0055731853) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.30243057) q[1];
sx q[1];
rz(-1.7453472) q[1];
sx q[1];
rz(-1.3329092) q[1];
rz(-pi) q[2];
x q[2];
rz(0.69435223) q[3];
sx q[3];
rz(-2.4789841) q[3];
sx q[3];
rz(-2.724805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.906189) q[2];
sx q[2];
rz(-1.0455422) q[2];
sx q[2];
rz(1.7380627) q[2];
rz(0.67052001) q[3];
sx q[3];
rz(-1.5454005) q[3];
sx q[3];
rz(-1.3129354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32187605) q[0];
sx q[0];
rz(-0.96915594) q[0];
sx q[0];
rz(-0.72625351) q[0];
rz(2.4655474) q[1];
sx q[1];
rz(-1.36422) q[1];
sx q[1];
rz(2.3647251) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1438734) q[0];
sx q[0];
rz(-1.6164427) q[0];
sx q[0];
rz(-0.93046988) q[0];
rz(2.9954722) q[2];
sx q[2];
rz(-1.6152799) q[2];
sx q[2];
rz(-1.5715949) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.54884262) q[1];
sx q[1];
rz(-1.934086) q[1];
sx q[1];
rz(2.8123358) q[1];
x q[2];
rz(-2.3915798) q[3];
sx q[3];
rz(-2.3976644) q[3];
sx q[3];
rz(2.9946526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1399416) q[2];
sx q[2];
rz(-0.17639128) q[2];
sx q[2];
rz(-1.4886935) q[2];
rz(2.0148924) q[3];
sx q[3];
rz(-1.9727581) q[3];
sx q[3];
rz(-0.53762976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5269932) q[0];
sx q[0];
rz(-2.499883) q[0];
sx q[0];
rz(-1.5747621) q[0];
rz(2.3552786) q[1];
sx q[1];
rz(-0.17335261) q[1];
sx q[1];
rz(-0.39549624) q[1];
rz(-2.4695071) q[2];
sx q[2];
rz(-2.0460049) q[2];
sx q[2];
rz(-2.6071215) q[2];
rz(2.366131) q[3];
sx q[3];
rz(-1.7453334) q[3];
sx q[3];
rz(-0.018669563) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
