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
rz(1.6504352) q[0];
sx q[0];
rz(2.7661134) q[0];
sx q[0];
rz(9.8162415) q[0];
rz(-2.904881) q[1];
sx q[1];
rz(-1.038329) q[1];
sx q[1];
rz(-0.84586039) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62097716) q[0];
sx q[0];
rz(-1.3987887) q[0];
sx q[0];
rz(0.67396236) q[0];
rz(-0.79997815) q[2];
sx q[2];
rz(-2.5350179) q[2];
sx q[2];
rz(2.6839566) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1497942) q[1];
sx q[1];
rz(-1.2956966) q[1];
sx q[1];
rz(0.28329785) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0835285) q[3];
sx q[3];
rz(-1.2635651) q[3];
sx q[3];
rz(1.787825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5819431) q[2];
sx q[2];
rz(-1.0407642) q[2];
sx q[2];
rz(-2.9259658) q[2];
rz(-0.98254472) q[3];
sx q[3];
rz(-1.6590786) q[3];
sx q[3];
rz(1.8854347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7243645) q[0];
sx q[0];
rz(-0.094014458) q[0];
sx q[0];
rz(-0.37658682) q[0];
rz(-1.6603893) q[1];
sx q[1];
rz(-1.9638289) q[1];
sx q[1];
rz(-2.1362163) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.004285) q[0];
sx q[0];
rz(-1.8958603) q[0];
sx q[0];
rz(-3.1336729) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2475904) q[2];
sx q[2];
rz(-1.6996222) q[2];
sx q[2];
rz(1.027838) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.59170106) q[1];
sx q[1];
rz(-1.5366067) q[1];
sx q[1];
rz(2.4728165) q[1];
x q[2];
rz(1.9485669) q[3];
sx q[3];
rz(-1.0527415) q[3];
sx q[3];
rz(-0.6686206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5873831) q[2];
sx q[2];
rz(-2.6770834) q[2];
sx q[2];
rz(-2.0655538) q[2];
rz(-0.46766034) q[3];
sx q[3];
rz(-0.83100072) q[3];
sx q[3];
rz(2.9369798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49863807) q[0];
sx q[0];
rz(-2.0641646) q[0];
sx q[0];
rz(-1.1109918) q[0];
rz(0.30771646) q[1];
sx q[1];
rz(-1.4765164) q[1];
sx q[1];
rz(-1.3002546) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0205524) q[0];
sx q[0];
rz(-0.55807923) q[0];
sx q[0];
rz(0.16772018) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5346393) q[2];
sx q[2];
rz(-2.1484408) q[2];
sx q[2];
rz(0.69427488) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.96836263) q[1];
sx q[1];
rz(-0.97665962) q[1];
sx q[1];
rz(2.6244375) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59107391) q[3];
sx q[3];
rz(-0.58286506) q[3];
sx q[3];
rz(1.7170563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51851455) q[2];
sx q[2];
rz(-1.0479505) q[2];
sx q[2];
rz(-2.9497362) q[2];
rz(2.0375371) q[3];
sx q[3];
rz(-0.42686978) q[3];
sx q[3];
rz(1.9328611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0722825) q[0];
sx q[0];
rz(-0.59232124) q[0];
sx q[0];
rz(0.88822547) q[0];
rz(-2.8690673) q[1];
sx q[1];
rz(-1.2963632) q[1];
sx q[1];
rz(1.9857508) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10242045) q[0];
sx q[0];
rz(-1.8646324) q[0];
sx q[0];
rz(0.0095902082) q[0];
rz(2.6510198) q[2];
sx q[2];
rz(-1.5340048) q[2];
sx q[2];
rz(1.8951777) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.72884293) q[1];
sx q[1];
rz(-0.83020617) q[1];
sx q[1];
rz(1.1918588) q[1];
rz(-0.4478332) q[3];
sx q[3];
rz(-1.574423) q[3];
sx q[3];
rz(1.8254395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3027975) q[2];
sx q[2];
rz(-2.8571627) q[2];
sx q[2];
rz(-2.3885041) q[2];
rz(-0.49790844) q[3];
sx q[3];
rz(-1.4830736) q[3];
sx q[3];
rz(-0.30043093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
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
rz(-0.65557757) q[0];
sx q[0];
rz(-2.1652554) q[0];
sx q[0];
rz(1.8560386) q[0];
rz(-1.9517508) q[1];
sx q[1];
rz(-2.4557476) q[1];
sx q[1];
rz(1.5815585) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.536895) q[0];
sx q[0];
rz(-1.492939) q[0];
sx q[0];
rz(-1.6998864) q[0];
rz(0.28376707) q[2];
sx q[2];
rz(-1.425425) q[2];
sx q[2];
rz(1.2710557) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.469857) q[1];
sx q[1];
rz(-1.3570045) q[1];
sx q[1];
rz(0.44624568) q[1];
x q[2];
rz(1.7722583) q[3];
sx q[3];
rz(-1.3115495) q[3];
sx q[3];
rz(-0.5539971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2994284) q[2];
sx q[2];
rz(-0.93141586) q[2];
sx q[2];
rz(3.0038339) q[2];
rz(2.6970741) q[3];
sx q[3];
rz(-1.7701745) q[3];
sx q[3];
rz(0.81454149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21823068) q[0];
sx q[0];
rz(-0.87319279) q[0];
sx q[0];
rz(-0.55106226) q[0];
rz(2.8463544) q[1];
sx q[1];
rz(-2.4199838) q[1];
sx q[1];
rz(2.3982184) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3819982) q[0];
sx q[0];
rz(-2.9659191) q[0];
sx q[0];
rz(0.86377527) q[0];
x q[1];
rz(-2.8001875) q[2];
sx q[2];
rz(-1.3479173) q[2];
sx q[2];
rz(-2.6472732) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36657143) q[1];
sx q[1];
rz(-1.4719392) q[1];
sx q[1];
rz(-0.84895433) q[1];
rz(-pi) q[2];
rz(0.94252276) q[3];
sx q[3];
rz(-0.97369558) q[3];
sx q[3];
rz(-1.3442849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.729852) q[2];
sx q[2];
rz(-2.9509632) q[2];
sx q[2];
rz(-2.8594678) q[2];
rz(2.7052346) q[3];
sx q[3];
rz(-1.8264419) q[3];
sx q[3];
rz(0.25562975) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9333164) q[0];
sx q[0];
rz(-1.0641119) q[0];
sx q[0];
rz(-2.5489885) q[0];
rz(-2.0049877) q[1];
sx q[1];
rz(-1.2717783) q[1];
sx q[1];
rz(-2.073854) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9737905) q[0];
sx q[0];
rz(-2.4616957) q[0];
sx q[0];
rz(0.25088422) q[0];
rz(2.8270193) q[2];
sx q[2];
rz(-2.7554263) q[2];
sx q[2];
rz(2.64606) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.065526) q[1];
sx q[1];
rz(-0.56030956) q[1];
sx q[1];
rz(2.1140772) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41530825) q[3];
sx q[3];
rz(-1.7303932) q[3];
sx q[3];
rz(0.94091614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.20256715) q[2];
sx q[2];
rz(-1.2747526) q[2];
sx q[2];
rz(2.1620046) q[2];
rz(0.12796399) q[3];
sx q[3];
rz(-0.94240677) q[3];
sx q[3];
rz(1.6894107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1487883) q[0];
sx q[0];
rz(-0.95071852) q[0];
sx q[0];
rz(0.08091452) q[0];
rz(0.12123904) q[1];
sx q[1];
rz(-0.67758766) q[1];
sx q[1];
rz(-2.0176719) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30643647) q[0];
sx q[0];
rz(-0.84266169) q[0];
sx q[0];
rz(-0.50424285) q[0];
rz(0.40718715) q[2];
sx q[2];
rz(-1.5296808) q[2];
sx q[2];
rz(-0.1830872) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2595265) q[1];
sx q[1];
rz(-1.9212544) q[1];
sx q[1];
rz(1.4653852) q[1];
rz(-0.1701351) q[3];
sx q[3];
rz(-1.0455344) q[3];
sx q[3];
rz(0.13290031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.27440444) q[2];
sx q[2];
rz(-0.97325456) q[2];
sx q[2];
rz(-1.493206) q[2];
rz(0.46227208) q[3];
sx q[3];
rz(-2.3267764) q[3];
sx q[3];
rz(1.4155242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82861519) q[0];
sx q[0];
rz(-1.7531489) q[0];
sx q[0];
rz(1.1653362) q[0];
rz(0.041821592) q[1];
sx q[1];
rz(-2.1526497) q[1];
sx q[1];
rz(1.8080541) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8671252) q[0];
sx q[0];
rz(-2.2750912) q[0];
sx q[0];
rz(-0.88732669) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.012537738) q[2];
sx q[2];
rz(-2.822753) q[2];
sx q[2];
rz(2.7818305) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0252777) q[1];
sx q[1];
rz(-1.1564558) q[1];
sx q[1];
rz(1.4715018) q[1];
rz(0.72527146) q[3];
sx q[3];
rz(-0.15003157) q[3];
sx q[3];
rz(-1.6440441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24451438) q[2];
sx q[2];
rz(-1.8316734) q[2];
sx q[2];
rz(-0.46373996) q[2];
rz(2.5044299) q[3];
sx q[3];
rz(-1.5183247) q[3];
sx q[3];
rz(2.7630828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8265182) q[0];
sx q[0];
rz(-1.0239064) q[0];
sx q[0];
rz(-1.6465323) q[0];
rz(-2.3118094) q[1];
sx q[1];
rz(-1.8779495) q[1];
sx q[1];
rz(-1.8216546) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49176106) q[0];
sx q[0];
rz(-1.6982838) q[0];
sx q[0];
rz(-2.7694635) q[0];
rz(-pi) q[1];
rz(2.7097335) q[2];
sx q[2];
rz(-1.3629902) q[2];
sx q[2];
rz(0.49655216) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4530427) q[1];
sx q[1];
rz(-0.31089766) q[1];
sx q[1];
rz(0.71605686) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0089031) q[3];
sx q[3];
rz(-1.3504538) q[3];
sx q[3];
rz(1.6957078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7497471) q[2];
sx q[2];
rz(-1.4186991) q[2];
sx q[2];
rz(0.5113655) q[2];
rz(-0.5558719) q[3];
sx q[3];
rz(-2.0248196) q[3];
sx q[3];
rz(-0.51978022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9075254) q[0];
sx q[0];
rz(-0.10744849) q[0];
sx q[0];
rz(0.43803122) q[0];
rz(-0.058649339) q[1];
sx q[1];
rz(-1.6390683) q[1];
sx q[1];
rz(-1.0674089) q[1];
rz(2.4406383) q[2];
sx q[2];
rz(-1.5264216) q[2];
sx q[2];
rz(3.0469372) q[2];
rz(-1.8773114) q[3];
sx q[3];
rz(-0.96592663) q[3];
sx q[3];
rz(2.929579) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
