OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4164299) q[0];
sx q[0];
rz(-0.13983146) q[0];
sx q[0];
rz(-2.5319985) q[0];
rz(0.66863376) q[1];
sx q[1];
rz(-2.2761087) q[1];
sx q[1];
rz(3.0545711) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4092769) q[0];
sx q[0];
rz(-1.9415932) q[0];
sx q[0];
rz(0.91863527) q[0];
rz(1.4423223) q[2];
sx q[2];
rz(-0.80768425) q[2];
sx q[2];
rz(-1.3070004) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.39057186) q[1];
sx q[1];
rz(-2.6387072) q[1];
sx q[1];
rz(1.8608171) q[1];
x q[2];
rz(3.0356785) q[3];
sx q[3];
rz(-2.7273791) q[3];
sx q[3];
rz(0.29990444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13880754) q[2];
sx q[2];
rz(-1.3802718) q[2];
sx q[2];
rz(0.37386093) q[2];
rz(-2.8047681) q[3];
sx q[3];
rz(-1.5461494) q[3];
sx q[3];
rz(-0.22836223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2574629) q[0];
sx q[0];
rz(-0.3733491) q[0];
sx q[0];
rz(1.194838) q[0];
rz(0.082611235) q[1];
sx q[1];
rz(-1.9742842) q[1];
sx q[1];
rz(-3.1412178) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99129358) q[0];
sx q[0];
rz(-0.55643686) q[0];
sx q[0];
rz(2.2668905) q[0];
x q[1];
rz(-1.3090918) q[2];
sx q[2];
rz(-1.3771025) q[2];
sx q[2];
rz(-1.2145834) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8723215) q[1];
sx q[1];
rz(-1.3355796) q[1];
sx q[1];
rz(-0.61551952) q[1];
rz(-pi) q[2];
rz(-0.11081103) q[3];
sx q[3];
rz(-0.46289819) q[3];
sx q[3];
rz(-3.0447931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.80883819) q[2];
sx q[2];
rz(-0.19503441) q[2];
sx q[2];
rz(-2.9648798) q[2];
rz(-0.79408944) q[3];
sx q[3];
rz(-2.3637171) q[3];
sx q[3];
rz(-0.55364048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2394543) q[0];
sx q[0];
rz(-2.4376526) q[0];
sx q[0];
rz(2.7368271) q[0];
rz(1.2813214) q[1];
sx q[1];
rz(-2.5679913) q[1];
sx q[1];
rz(1.8331029) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93861261) q[0];
sx q[0];
rz(-0.70883195) q[0];
sx q[0];
rz(2.3908486) q[0];
x q[1];
rz(-2.6042125) q[2];
sx q[2];
rz(-1.7747305) q[2];
sx q[2];
rz(0.98762074) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.96924671) q[1];
sx q[1];
rz(-1.3174651) q[1];
sx q[1];
rz(2.2971056) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6609459) q[3];
sx q[3];
rz(-2.0079068) q[3];
sx q[3];
rz(2.3428832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9222766) q[2];
sx q[2];
rz(-2.6185161) q[2];
sx q[2];
rz(3.1211839) q[2];
rz(1.071788) q[3];
sx q[3];
rz(-1.1276378) q[3];
sx q[3];
rz(-0.66334692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14169176) q[0];
sx q[0];
rz(-0.90943709) q[0];
sx q[0];
rz(0.91598696) q[0];
rz(0.46332106) q[1];
sx q[1];
rz(-2.0756192) q[1];
sx q[1];
rz(-1.0571009) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0650478) q[0];
sx q[0];
rz(-1.1335982) q[0];
sx q[0];
rz(-1.5426427) q[0];
x q[1];
rz(-2.5713423) q[2];
sx q[2];
rz(-2.1006561) q[2];
sx q[2];
rz(-0.8023163) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.033356655) q[1];
sx q[1];
rz(-1.557096) q[1];
sx q[1];
rz(-1.800888) q[1];
rz(1.3435059) q[3];
sx q[3];
rz(-0.82633457) q[3];
sx q[3];
rz(2.8518761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.89381924) q[2];
sx q[2];
rz(-0.63988581) q[2];
sx q[2];
rz(0.34269732) q[2];
rz(-1.6977067) q[3];
sx q[3];
rz(-2.2659437) q[3];
sx q[3];
rz(-2.4263583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7856359) q[0];
sx q[0];
rz(-0.56348339) q[0];
sx q[0];
rz(-3.0551531) q[0];
rz(-1.7516288) q[1];
sx q[1];
rz(-0.93170634) q[1];
sx q[1];
rz(-0.46868971) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8353032) q[0];
sx q[0];
rz(-1.1790743) q[0];
sx q[0];
rz(-0.83175559) q[0];
x q[1];
rz(-0.52207559) q[2];
sx q[2];
rz(-2.7895088) q[2];
sx q[2];
rz(3.0662231) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9710755) q[1];
sx q[1];
rz(-1.7726388) q[1];
sx q[1];
rz(1.7032196) q[1];
rz(-pi) q[2];
rz(1.4851941) q[3];
sx q[3];
rz(-0.75827956) q[3];
sx q[3];
rz(-3.0183834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.1420574) q[2];
sx q[2];
rz(-0.88025981) q[2];
sx q[2];
rz(-2.2772677) q[2];
rz(-2.9243829) q[3];
sx q[3];
rz(-0.61621284) q[3];
sx q[3];
rz(-0.29278452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42751673) q[0];
sx q[0];
rz(-1.7345411) q[0];
sx q[0];
rz(2.9445904) q[0];
rz(-1.3621832) q[1];
sx q[1];
rz(-1.4809337) q[1];
sx q[1];
rz(-0.33624712) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7173548) q[0];
sx q[0];
rz(-0.80168085) q[0];
sx q[0];
rz(2.5108811) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3669846) q[2];
sx q[2];
rz(-2.3850072) q[2];
sx q[2];
rz(-0.87385439) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.058192249) q[1];
sx q[1];
rz(-2.5501049) q[1];
sx q[1];
rz(3.0565492) q[1];
rz(-pi) q[2];
rz(-1.398596) q[3];
sx q[3];
rz(-2.2599054) q[3];
sx q[3];
rz(-2.8840051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6062935) q[2];
sx q[2];
rz(-1.4761816) q[2];
sx q[2];
rz(2.2952648) q[2];
rz(1.2396631) q[3];
sx q[3];
rz(-1.2318434) q[3];
sx q[3];
rz(-3.1404176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7476615) q[0];
sx q[0];
rz(-2.1037536) q[0];
sx q[0];
rz(-2.8826662) q[0];
rz(1.3461643) q[1];
sx q[1];
rz(-1.7566453) q[1];
sx q[1];
rz(-1.0940201) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0897113) q[0];
sx q[0];
rz(-1.9117172) q[0];
sx q[0];
rz(1.107035) q[0];
rz(-0.012374087) q[2];
sx q[2];
rz(-1.7933328) q[2];
sx q[2];
rz(-0.48977938) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3257521) q[1];
sx q[1];
rz(-0.098512352) q[1];
sx q[1];
rz(-1.8580395) q[1];
x q[2];
rz(2.9225227) q[3];
sx q[3];
rz(-3.0411358) q[3];
sx q[3];
rz(0.62517525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.1743494) q[2];
sx q[2];
rz(-1.3568342) q[2];
sx q[2];
rz(-2.8209177) q[2];
rz(2.705412) q[3];
sx q[3];
rz(-2.6685721) q[3];
sx q[3];
rz(2.6935553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2748579) q[0];
sx q[0];
rz(-0.47645706) q[0];
sx q[0];
rz(2.136769) q[0];
rz(2.5514305) q[1];
sx q[1];
rz(-2.2361123) q[1];
sx q[1];
rz(2.337713) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.852254) q[0];
sx q[0];
rz(-0.58710557) q[0];
sx q[0];
rz(3.0703074) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5875823) q[2];
sx q[2];
rz(-1.4129352) q[2];
sx q[2];
rz(2.6487034) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.67919532) q[1];
sx q[1];
rz(-3.052366) q[1];
sx q[1];
rz(1.4560844) q[1];
rz(2.478066) q[3];
sx q[3];
rz(-1.4995721) q[3];
sx q[3];
rz(1.0458667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.33264318) q[2];
sx q[2];
rz(-0.84471622) q[2];
sx q[2];
rz(2.1311029) q[2];
rz(-0.60339749) q[3];
sx q[3];
rz(-1.5248652) q[3];
sx q[3];
rz(0.55019125) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6036966) q[0];
sx q[0];
rz(-1.864707) q[0];
sx q[0];
rz(-0.4075152) q[0];
rz(-2.852476) q[1];
sx q[1];
rz(-1.1228077) q[1];
sx q[1];
rz(2.3908652) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6347046) q[0];
sx q[0];
rz(-1.2872211) q[0];
sx q[0];
rz(-0.30446913) q[0];
rz(1.4892011) q[2];
sx q[2];
rz(-1.8494693) q[2];
sx q[2];
rz(-0.96688731) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6299906) q[1];
sx q[1];
rz(-1.8774319) q[1];
sx q[1];
rz(-2.929045) q[1];
rz(-pi) q[2];
rz(3.0575033) q[3];
sx q[3];
rz(-0.30232271) q[3];
sx q[3];
rz(-1.7526527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0728545) q[2];
sx q[2];
rz(-0.72454238) q[2];
sx q[2];
rz(-1.9753974) q[2];
rz(1.3646305) q[3];
sx q[3];
rz(-2.3659673) q[3];
sx q[3];
rz(-0.021818074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5584548) q[0];
sx q[0];
rz(-2.3151509) q[0];
sx q[0];
rz(-1.3903842) q[0];
rz(-0.33069262) q[1];
sx q[1];
rz(-2.3762517) q[1];
sx q[1];
rz(-1.4601382) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6351785) q[0];
sx q[0];
rz(-1.9902475) q[0];
sx q[0];
rz(1.0181396) q[0];
rz(-pi) q[1];
rz(0.023256217) q[2];
sx q[2];
rz(-2.0746982) q[2];
sx q[2];
rz(1.1018673) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8148988) q[1];
sx q[1];
rz(-1.1732475) q[1];
sx q[1];
rz(-1.7705998) q[1];
rz(-pi) q[2];
rz(2.2721685) q[3];
sx q[3];
rz(-2.2079909) q[3];
sx q[3];
rz(-0.77995342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7918487) q[2];
sx q[2];
rz(-1.3246374) q[2];
sx q[2];
rz(-1.4617408) q[2];
rz(-2.0215624) q[3];
sx q[3];
rz(-0.62168613) q[3];
sx q[3];
rz(-2.6859443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.6256975) q[0];
sx q[0];
rz(-1.5213756) q[0];
sx q[0];
rz(-1.9949927) q[0];
rz(1.3810146) q[1];
sx q[1];
rz(-1.2881423) q[1];
sx q[1];
rz(-1.2013411) q[1];
rz(-2.7881651) q[2];
sx q[2];
rz(-2.3218487) q[2];
sx q[2];
rz(-0.29816366) q[2];
rz(1.9588884) q[3];
sx q[3];
rz(-0.93508616) q[3];
sx q[3];
rz(-1.3133776) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
