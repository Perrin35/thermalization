OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.779939) q[0];
sx q[0];
rz(-0.10804636) q[0];
sx q[0];
rz(1.894423) q[0];
rz(-0.78020686) q[1];
sx q[1];
rz(-1.127004) q[1];
sx q[1];
rz(-0.74581528) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5818258) q[0];
sx q[0];
rz(-0.39224658) q[0];
sx q[0];
rz(2.0779209) q[0];
rz(-0.0086486625) q[2];
sx q[2];
rz(-1.7315355) q[2];
sx q[2];
rz(2.6313105) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8565897) q[1];
sx q[1];
rz(-1.7816646) q[1];
sx q[1];
rz(3.0233356) q[1];
x q[2];
rz(2.8549544) q[3];
sx q[3];
rz(-2.6034546) q[3];
sx q[3];
rz(-1.2179483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1625533) q[2];
sx q[2];
rz(-2.5472842) q[2];
sx q[2];
rz(1.7621367) q[2];
rz(0.72921324) q[3];
sx q[3];
rz(-2.3766434) q[3];
sx q[3];
rz(-2.2123857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.5586229) q[0];
sx q[0];
rz(-2.0780777) q[0];
sx q[0];
rz(0.24120086) q[0];
rz(0.25602117) q[1];
sx q[1];
rz(-1.0216917) q[1];
sx q[1];
rz(2.3604438) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1560416) q[0];
sx q[0];
rz(-2.3325778) q[0];
sx q[0];
rz(2.3963905) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5767869) q[2];
sx q[2];
rz(-1.6132659) q[2];
sx q[2];
rz(-3.0926585) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.95805146) q[1];
sx q[1];
rz(-2.1528917) q[1];
sx q[1];
rz(0.84611012) q[1];
rz(1.2595176) q[3];
sx q[3];
rz(-2.9972509) q[3];
sx q[3];
rz(0.33931574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8476734) q[2];
sx q[2];
rz(-1.5896229) q[2];
sx q[2];
rz(1.5839362) q[2];
rz(-3.1368351) q[3];
sx q[3];
rz(-0.59499756) q[3];
sx q[3];
rz(2.7958272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7541589) q[0];
sx q[0];
rz(-1.6922981) q[0];
sx q[0];
rz(-2.9505728) q[0];
rz(-0.35107958) q[1];
sx q[1];
rz(-0.43126884) q[1];
sx q[1];
rz(-1.692159) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3905778) q[0];
sx q[0];
rz(-1.6216941) q[0];
sx q[0];
rz(-0.33945947) q[0];
x q[1];
rz(0.44523737) q[2];
sx q[2];
rz(-1.9539333) q[2];
sx q[2];
rz(-0.68510011) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8061132) q[1];
sx q[1];
rz(-1.8361143) q[1];
sx q[1];
rz(0.47228864) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5496629) q[3];
sx q[3];
rz(-1.283965) q[3];
sx q[3];
rz(-1.5805707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8423975) q[2];
sx q[2];
rz(-1.4292052) q[2];
sx q[2];
rz(3.0136285) q[2];
rz(1.1658824) q[3];
sx q[3];
rz(-0.88763014) q[3];
sx q[3];
rz(0.33795801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8562451) q[0];
sx q[0];
rz(-2.0576394) q[0];
sx q[0];
rz(1.4803084) q[0];
rz(3.0604494) q[1];
sx q[1];
rz(-2.0442043) q[1];
sx q[1];
rz(-2.2611484) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1937516) q[0];
sx q[0];
rz(-1.546085) q[0];
sx q[0];
rz(-1.5642883) q[0];
x q[1];
rz(-1.501785) q[2];
sx q[2];
rz(-2.5216148) q[2];
sx q[2];
rz(2.2295956) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.26759155) q[1];
sx q[1];
rz(-2.1028215) q[1];
sx q[1];
rz(1.6997972) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.038136884) q[3];
sx q[3];
rz(-2.5835369) q[3];
sx q[3];
rz(1.3636774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.60260281) q[2];
sx q[2];
rz(-0.67402855) q[2];
sx q[2];
rz(-1.7146141) q[2];
rz(-1.9849298) q[3];
sx q[3];
rz(-1.7159228) q[3];
sx q[3];
rz(-0.15779933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91931525) q[0];
sx q[0];
rz(-0.8150402) q[0];
sx q[0];
rz(-2.2015233) q[0];
rz(0.4982416) q[1];
sx q[1];
rz(-0.91612852) q[1];
sx q[1];
rz(-2.0169651) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2021411) q[0];
sx q[0];
rz(-1.5312342) q[0];
sx q[0];
rz(0.8673773) q[0];
rz(-pi) q[1];
rz(1.1246936) q[2];
sx q[2];
rz(-0.44844018) q[2];
sx q[2];
rz(-0.61277991) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4044721) q[1];
sx q[1];
rz(-1.1386765) q[1];
sx q[1];
rz(-2.4678556) q[1];
rz(-pi) q[2];
rz(0.78557555) q[3];
sx q[3];
rz(-1.5858558) q[3];
sx q[3];
rz(2.756898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4336832) q[2];
sx q[2];
rz(-2.3072672) q[2];
sx q[2];
rz(1.4617807) q[2];
rz(2.8953569) q[3];
sx q[3];
rz(-0.88053954) q[3];
sx q[3];
rz(-1.0730526) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4435302) q[0];
sx q[0];
rz(-1.9140697) q[0];
sx q[0];
rz(0.45502934) q[0];
rz(-2.9235234) q[1];
sx q[1];
rz(-2.2360305) q[1];
sx q[1];
rz(0.47854447) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8826101) q[0];
sx q[0];
rz(-2.4702669) q[0];
sx q[0];
rz(0.15780003) q[0];
rz(-pi) q[1];
rz(2.6784513) q[2];
sx q[2];
rz(-0.69370334) q[2];
sx q[2];
rz(2.5460473) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8560972) q[1];
sx q[1];
rz(-1.2427274) q[1];
sx q[1];
rz(2.0542007) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4535937) q[3];
sx q[3];
rz(-1.08537) q[3];
sx q[3];
rz(-1.7958876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9559418) q[2];
sx q[2];
rz(-1.4943244) q[2];
sx q[2];
rz(0.68598023) q[2];
rz(1.8020804) q[3];
sx q[3];
rz(-1.1687665) q[3];
sx q[3];
rz(-1.7495988) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5865536) q[0];
sx q[0];
rz(-1.3864484) q[0];
sx q[0];
rz(-2.6626124) q[0];
rz(2.2827177) q[1];
sx q[1];
rz(-2.5996467) q[1];
sx q[1];
rz(-0.10471334) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0901439) q[0];
sx q[0];
rz(-1.3810754) q[0];
sx q[0];
rz(-0.39162292) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62868406) q[2];
sx q[2];
rz(-2.5631944) q[2];
sx q[2];
rz(-2.4378928) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5471902) q[1];
sx q[1];
rz(-2.5195751) q[1];
sx q[1];
rz(1.1941431) q[1];
rz(2.0488304) q[3];
sx q[3];
rz(-1.6217188) q[3];
sx q[3];
rz(-1.6855406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0571478) q[2];
sx q[2];
rz(-1.2556602) q[2];
sx q[2];
rz(2.2596333) q[2];
rz(3.0706578) q[3];
sx q[3];
rz(-1.2627914) q[3];
sx q[3];
rz(-0.96496636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98081812) q[0];
sx q[0];
rz(-2.6698298) q[0];
sx q[0];
rz(-1.2534575) q[0];
rz(-2.4313633) q[1];
sx q[1];
rz(-1.9435147) q[1];
sx q[1];
rz(2.0517147) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80515758) q[0];
sx q[0];
rz(-1.1602279) q[0];
sx q[0];
rz(0.76469501) q[0];
rz(-pi) q[1];
rz(-2.5023191) q[2];
sx q[2];
rz(-2.6007923) q[2];
sx q[2];
rz(1.6317473) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8556887) q[1];
sx q[1];
rz(-0.82997417) q[1];
sx q[1];
rz(-0.94372933) q[1];
x q[2];
rz(-2.9097611) q[3];
sx q[3];
rz(-2.5527195) q[3];
sx q[3];
rz(0.26118054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.39500427) q[2];
sx q[2];
rz(-2.2825664) q[2];
sx q[2];
rz(1.0224226) q[2];
rz(2.5452781) q[3];
sx q[3];
rz(-2.3583581) q[3];
sx q[3];
rz(0.35308009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31350598) q[0];
sx q[0];
rz(-2.6441898) q[0];
sx q[0];
rz(-1.5561546) q[0];
rz(1.4900788) q[1];
sx q[1];
rz(-0.45526344) q[1];
sx q[1];
rz(-2.1447287) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6108902) q[0];
sx q[0];
rz(-0.35323745) q[0];
sx q[0];
rz(-0.20521407) q[0];
rz(-1.385594) q[2];
sx q[2];
rz(-1.569199) q[2];
sx q[2];
rz(-1.314437) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.026406757) q[1];
sx q[1];
rz(-0.86167012) q[1];
sx q[1];
rz(0.89964189) q[1];
rz(-pi) q[2];
rz(2.6410651) q[3];
sx q[3];
rz(-1.0382892) q[3];
sx q[3];
rz(2.4103269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8263714) q[2];
sx q[2];
rz(-2.4085277) q[2];
sx q[2];
rz(-0.23615393) q[2];
rz(0.30596966) q[3];
sx q[3];
rz(-1.1924084) q[3];
sx q[3];
rz(1.6570305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21208256) q[0];
sx q[0];
rz(-0.65123737) q[0];
sx q[0];
rz(-0.65810743) q[0];
rz(1.8923538) q[1];
sx q[1];
rz(-0.58964261) q[1];
sx q[1];
rz(0.49159893) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99557331) q[0];
sx q[0];
rz(-0.73500145) q[0];
sx q[0];
rz(-2.4400178) q[0];
x q[1];
rz(1.3479606) q[2];
sx q[2];
rz(-1.0165983) q[2];
sx q[2];
rz(1.8775307) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.39855865) q[1];
sx q[1];
rz(-1.5614206) q[1];
sx q[1];
rz(-0.27573632) q[1];
rz(-pi) q[2];
rz(-1.0820457) q[3];
sx q[3];
rz(-1.6433892) q[3];
sx q[3];
rz(-1.5047764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3601941) q[2];
sx q[2];
rz(-2.7935544) q[2];
sx q[2];
rz(-1.4813102) q[2];
rz(2.4252452) q[3];
sx q[3];
rz(-2.4951388) q[3];
sx q[3];
rz(-0.20814482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3502055) q[0];
sx q[0];
rz(-1.5443784) q[0];
sx q[0];
rz(-2.7271893) q[0];
rz(2.1898337) q[1];
sx q[1];
rz(-1.2928243) q[1];
sx q[1];
rz(-1.181319) q[1];
rz(-0.70474456) q[2];
sx q[2];
rz(-0.8316883) q[2];
sx q[2];
rz(-1.4744454) q[2];
rz(-2.9708859) q[3];
sx q[3];
rz(-2.2424663) q[3];
sx q[3];
rz(-0.57244931) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
