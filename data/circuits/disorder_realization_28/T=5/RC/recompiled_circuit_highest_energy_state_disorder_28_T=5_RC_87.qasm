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
rz(-2.7653341) q[1];
sx q[1];
rz(-0.19317746) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8602919) q[0];
sx q[0];
rz(-2.0112462) q[0];
sx q[0];
rz(2.4971636) q[0];
rz(-pi) q[1];
rz(-0.12234064) q[2];
sx q[2];
rz(-2.2723778) q[2];
sx q[2];
rz(2.0035494) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.2235893) q[1];
rz(-1.3985996) q[3];
sx q[3];
rz(-1.2463257) q[3];
sx q[3];
rz(-2.8419122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2304077) q[2];
sx q[2];
rz(-2.513803) q[2];
sx q[2];
rz(0.5578624) q[2];
rz(1.9134391) q[3];
sx q[3];
rz(-1.4598673) q[3];
sx q[3];
rz(0.094999878) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2896344) q[0];
sx q[0];
rz(-1.292922) q[0];
sx q[0];
rz(-1.5767545) q[0];
rz(1.9948888) q[1];
sx q[1];
rz(-1.7161938) q[1];
sx q[1];
rz(1.0316521) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0609425) q[0];
sx q[0];
rz(-0.58851465) q[0];
sx q[0];
rz(1.1400162) q[0];
rz(-0.41298683) q[2];
sx q[2];
rz(-1.7888165) q[2];
sx q[2];
rz(3.060201) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.37687846) q[1];
sx q[1];
rz(-2.2441314) q[1];
sx q[1];
rz(-0.44133472) q[1];
rz(-pi) q[2];
rz(-0.50526039) q[3];
sx q[3];
rz(-1.7499515) q[3];
sx q[3];
rz(-3.1324838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8138294) q[2];
sx q[2];
rz(-2.8507865) q[2];
sx q[2];
rz(1.4168868) q[2];
rz(-2.318577) q[3];
sx q[3];
rz(-1.5282642) q[3];
sx q[3];
rz(2.4970162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7773892) q[0];
sx q[0];
rz(-1.1193898) q[0];
sx q[0];
rz(2.7893344) q[0];
rz(2.7525355) q[1];
sx q[1];
rz(-1.9720826) q[1];
sx q[1];
rz(0.22638098) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8733559) q[0];
sx q[0];
rz(-1.4882636) q[0];
sx q[0];
rz(-1.8186452) q[0];
rz(-2.1445072) q[2];
sx q[2];
rz(-1.2752394) q[2];
sx q[2];
rz(-1.3546582) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6701874) q[1];
sx q[1];
rz(-2.4627663) q[1];
sx q[1];
rz(0.77969867) q[1];
x q[2];
rz(1.6466633) q[3];
sx q[3];
rz(-1.1359204) q[3];
sx q[3];
rz(-2.7135389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9809197) q[2];
sx q[2];
rz(-1.0272762) q[2];
sx q[2];
rz(0.36160198) q[2];
rz(2.8412039) q[3];
sx q[3];
rz(-1.8301679) q[3];
sx q[3];
rz(-0.96295199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2196197) q[0];
sx q[0];
rz(-2.0491056) q[0];
sx q[0];
rz(3.0203982) q[0];
rz(1.1990064) q[1];
sx q[1];
rz(-1.0795178) q[1];
sx q[1];
rz(-1.2746864) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.072678) q[0];
sx q[0];
rz(-1.5240178) q[0];
sx q[0];
rz(-0.22225265) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.936144) q[2];
sx q[2];
rz(-1.5549855) q[2];
sx q[2];
rz(-0.39218536) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1772985) q[1];
sx q[1];
rz(-1.3508479) q[1];
sx q[1];
rz(-2.5581058) q[1];
rz(-pi) q[2];
rz(-1.9363112) q[3];
sx q[3];
rz(-1.707122) q[3];
sx q[3];
rz(2.8773235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.092209665) q[2];
sx q[2];
rz(-1.4468687) q[2];
sx q[2];
rz(1.764074) q[2];
rz(-1.1285909) q[3];
sx q[3];
rz(-2.1994574) q[3];
sx q[3];
rz(-0.82730627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89151299) q[0];
sx q[0];
rz(-2.2703607) q[0];
sx q[0];
rz(-1.0846035) q[0];
rz(-2.4705823) q[1];
sx q[1];
rz(-1.6295461) q[1];
sx q[1];
rz(3.0674518) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29956578) q[0];
sx q[0];
rz(-1.8447829) q[0];
sx q[0];
rz(-0.58273231) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2535278) q[2];
sx q[2];
rz(-1.2937577) q[2];
sx q[2];
rz(-0.46169567) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.45696631) q[1];
sx q[1];
rz(-1.8408482) q[1];
sx q[1];
rz(-0.96182386) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0054194) q[3];
sx q[3];
rz(-1.2981997) q[3];
sx q[3];
rz(-2.5933602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.025617754) q[2];
sx q[2];
rz(-2.3927549) q[2];
sx q[2];
rz(-0.96308127) q[2];
rz(1.3747619) q[3];
sx q[3];
rz(-1.129496) q[3];
sx q[3];
rz(-2.6095384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2087723) q[0];
sx q[0];
rz(-2.8498579) q[0];
sx q[0];
rz(-1.2579086) q[0];
rz(-2.3014297) q[1];
sx q[1];
rz(-1.9636619) q[1];
sx q[1];
rz(-2.449583) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7888579) q[0];
sx q[0];
rz(-1.8957355) q[0];
sx q[0];
rz(-2.8558225) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18628405) q[2];
sx q[2];
rz(-2.3669503) q[2];
sx q[2];
rz(1.2828522) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4547448) q[1];
sx q[1];
rz(-0.94986373) q[1];
sx q[1];
rz(-2.5140155) q[1];
rz(2.5049772) q[3];
sx q[3];
rz(-1.1061543) q[3];
sx q[3];
rz(-0.98004228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.96364) q[2];
sx q[2];
rz(-1.0963564) q[2];
sx q[2];
rz(0.5160416) q[2];
rz(1.3456723) q[3];
sx q[3];
rz(-2.7005152) q[3];
sx q[3];
rz(-1.0015944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26821414) q[0];
sx q[0];
rz(-2.8554947) q[0];
sx q[0];
rz(2.1116665) q[0];
rz(2.4870807) q[1];
sx q[1];
rz(-1.8067358) q[1];
sx q[1];
rz(-2.9676504) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2692341) q[0];
sx q[0];
rz(-2.1790811) q[0];
sx q[0];
rz(-2.48003) q[0];
rz(-0.80751626) q[2];
sx q[2];
rz(-1.9587729) q[2];
sx q[2];
rz(-2.9569851) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0566114) q[1];
sx q[1];
rz(-2.1826996) q[1];
sx q[1];
rz(-0.33501825) q[1];
x q[2];
rz(1.8344363) q[3];
sx q[3];
rz(-1.0705433) q[3];
sx q[3];
rz(-1.6135581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.30764636) q[2];
sx q[2];
rz(-0.53457326) q[2];
sx q[2];
rz(1.3391116) q[2];
rz(2.3329959) q[3];
sx q[3];
rz(-1.1491038) q[3];
sx q[3];
rz(-1.8554525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22063743) q[0];
sx q[0];
rz(-2.0957102) q[0];
sx q[0];
rz(-1.9761696) q[0];
rz(0.19365817) q[1];
sx q[1];
rz(-0.8332738) q[1];
sx q[1];
rz(-1.1509034) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5315463) q[0];
sx q[0];
rz(-2.4360516) q[0];
sx q[0];
rz(0.39682589) q[0];
rz(-pi) q[1];
rz(-0.1702383) q[2];
sx q[2];
rz(-0.81184978) q[2];
sx q[2];
rz(-1.591631) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2886758) q[1];
sx q[1];
rz(-1.9948469) q[1];
sx q[1];
rz(2.1217594) q[1];
rz(-pi) q[2];
rz(2.1459313) q[3];
sx q[3];
rz(-1.0049977) q[3];
sx q[3];
rz(3.107323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99106818) q[0];
sx q[0];
rz(-3.0758698) q[0];
sx q[0];
rz(-2.7142628) q[0];
rz(1.4409675) q[1];
sx q[1];
rz(-1.7040375) q[1];
sx q[1];
rz(0.89868054) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3220427) q[0];
sx q[0];
rz(-1.2614095) q[0];
sx q[0];
rz(2.0218819) q[0];
x q[1];
rz(-0.9814453) q[2];
sx q[2];
rz(-1.169983) q[2];
sx q[2];
rz(-3.1360195) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2262804) q[1];
sx q[1];
rz(-1.8050005) q[1];
sx q[1];
rz(-2.9620902) q[1];
x q[2];
rz(-2.4472404) q[3];
sx q[3];
rz(-0.66260856) q[3];
sx q[3];
rz(2.724805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.906189) q[2];
sx q[2];
rz(-1.0455422) q[2];
sx q[2];
rz(1.7380627) q[2];
rz(-0.67052001) q[3];
sx q[3];
rz(-1.5961921) q[3];
sx q[3];
rz(1.8286573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.7773726) q[1];
sx q[1];
rz(-2.3647251) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51190864) q[0];
sx q[0];
rz(-0.64172318) q[0];
sx q[0];
rz(-1.4944906) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6157584) q[2];
sx q[2];
rz(-1.7167712) q[2];
sx q[2];
rz(-0.0073429664) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.59275) q[1];
sx q[1];
rz(-1.934086) q[1];
sx q[1];
rz(2.8123358) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75001287) q[3];
sx q[3];
rz(-2.3976644) q[3];
sx q[3];
rz(2.9946526) q[3];
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
rz(-1.9727581) q[3];
sx q[3];
rz(2.6039629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61459944) q[0];
sx q[0];
rz(-0.64170964) q[0];
sx q[0];
rz(1.5668305) q[0];
rz(0.78631403) q[1];
sx q[1];
rz(-2.96824) q[1];
sx q[1];
rz(2.7460964) q[1];
rz(0.98914374) q[2];
sx q[2];
rz(-2.157515) q[2];
sx q[2];
rz(-0.6872057) q[2];
rz(-1.3287104) q[3];
sx q[3];
rz(-2.3314706) q[3];
sx q[3];
rz(-1.7580845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
