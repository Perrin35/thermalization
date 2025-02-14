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
rz(-1.2404233) q[0];
sx q[0];
rz(-1.197553) q[0];
sx q[0];
rz(-0.21633202) q[0];
rz(0.95739111) q[1];
sx q[1];
rz(-2.4654145) q[1];
sx q[1];
rz(-1.6300936) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0656242) q[0];
sx q[0];
rz(-2.127021) q[0];
sx q[0];
rz(-1.5821304) q[0];
x q[1];
rz(-1.2443107) q[2];
sx q[2];
rz(-2.0225581) q[2];
sx q[2];
rz(1.8725852) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.18876702) q[1];
sx q[1];
rz(-0.61453846) q[1];
sx q[1];
rz(-0.60093083) q[1];
rz(-1.1153631) q[3];
sx q[3];
rz(-2.3854227) q[3];
sx q[3];
rz(-0.84634534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1067074) q[2];
sx q[2];
rz(-0.35197508) q[2];
sx q[2];
rz(2.0931639) q[2];
rz(0.18167051) q[3];
sx q[3];
rz(-2.1766267) q[3];
sx q[3];
rz(2.488193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.492391) q[0];
sx q[0];
rz(-0.94434706) q[0];
sx q[0];
rz(0.44678584) q[0];
rz(-0.88042879) q[1];
sx q[1];
rz(-1.7767521) q[1];
sx q[1];
rz(0.78278881) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7004958) q[0];
sx q[0];
rz(-1.5095599) q[0];
sx q[0];
rz(-1.3260617) q[0];
rz(-0.57116429) q[2];
sx q[2];
rz(-2.4413902) q[2];
sx q[2];
rz(-0.68298662) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.70322733) q[1];
sx q[1];
rz(-1.2709054) q[1];
sx q[1];
rz(-0.48734003) q[1];
rz(-pi) q[2];
rz(2.2530678) q[3];
sx q[3];
rz(-1.9491213) q[3];
sx q[3];
rz(2.494787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.030152628) q[2];
sx q[2];
rz(-1.4806662) q[2];
sx q[2];
rz(0.97935575) q[2];
rz(2.7566946) q[3];
sx q[3];
rz(-1.220547) q[3];
sx q[3];
rz(0.16429193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48591831) q[0];
sx q[0];
rz(-3.0874708) q[0];
sx q[0];
rz(2.3517877) q[0];
rz(2.9549331) q[1];
sx q[1];
rz(-1.7402382) q[1];
sx q[1];
rz(0.99367118) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6826166) q[0];
sx q[0];
rz(-1.0280711) q[0];
sx q[0];
rz(-1.8629462) q[0];
x q[1];
rz(-1.7970071) q[2];
sx q[2];
rz(-2.1717584) q[2];
sx q[2];
rz(0.41659875) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4504272) q[1];
sx q[1];
rz(-1.1466999) q[1];
sx q[1];
rz(2.0753839) q[1];
rz(0.72088269) q[3];
sx q[3];
rz(-1.7768806) q[3];
sx q[3];
rz(2.4997847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8197202) q[2];
sx q[2];
rz(-1.1178144) q[2];
sx q[2];
rz(2.7703088) q[2];
rz(-2.7927981) q[3];
sx q[3];
rz(-1.0943509) q[3];
sx q[3];
rz(0.5955407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.45469859) q[0];
sx q[0];
rz(-1.0201539) q[0];
sx q[0];
rz(1.7373079) q[0];
rz(-2.4644201) q[1];
sx q[1];
rz(-1.9858457) q[1];
sx q[1];
rz(1.6489395) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20805222) q[0];
sx q[0];
rz(-1.3553936) q[0];
sx q[0];
rz(-2.893704) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1588232) q[2];
sx q[2];
rz(-2.6361536) q[2];
sx q[2];
rz(0.44014058) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.75912913) q[1];
sx q[1];
rz(-1.3037523) q[1];
sx q[1];
rz(-1.8014531) q[1];
rz(-2.2954313) q[3];
sx q[3];
rz(-1.3337047) q[3];
sx q[3];
rz(2.2711636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6877785) q[2];
sx q[2];
rz(-2.1731589) q[2];
sx q[2];
rz(-0.34823927) q[2];
rz(1.4549152) q[3];
sx q[3];
rz(-1.429052) q[3];
sx q[3];
rz(1.3454364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1763879) q[0];
sx q[0];
rz(-2.5768953) q[0];
sx q[0];
rz(0.89214605) q[0];
rz(-2.6773894) q[1];
sx q[1];
rz(-1.8966388) q[1];
sx q[1];
rz(-1.8468599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87128042) q[0];
sx q[0];
rz(-2.0593606) q[0];
sx q[0];
rz(2.110092) q[0];
x q[1];
rz(1.4424397) q[2];
sx q[2];
rz(-2.1048628) q[2];
sx q[2];
rz(-0.72557025) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0820992) q[1];
sx q[1];
rz(-2.0096742) q[1];
sx q[1];
rz(2.782269) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2512359) q[3];
sx q[3];
rz(-2.0687752) q[3];
sx q[3];
rz(2.339956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3225473) q[2];
sx q[2];
rz(-2.5307541) q[2];
sx q[2];
rz(-2.5757705) q[2];
rz(3.0602509) q[3];
sx q[3];
rz(-2.1606725) q[3];
sx q[3];
rz(-2.5210023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.464798) q[0];
sx q[0];
rz(-0.066635266) q[0];
sx q[0];
rz(-1.5555405) q[0];
rz(1.0606891) q[1];
sx q[1];
rz(-1.5763177) q[1];
sx q[1];
rz(0.63180822) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6752967) q[0];
sx q[0];
rz(-1.491437) q[0];
sx q[0];
rz(-0.23510374) q[0];
rz(-1.00497) q[2];
sx q[2];
rz(-2.3477051) q[2];
sx q[2];
rz(-2.6086406) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0630707) q[1];
sx q[1];
rz(-1.0390738) q[1];
sx q[1];
rz(-2.7460981) q[1];
rz(-pi) q[2];
rz(-1.4238213) q[3];
sx q[3];
rz(-0.89255652) q[3];
sx q[3];
rz(-3.105046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1598728) q[2];
sx q[2];
rz(-2.0231415) q[2];
sx q[2];
rz(-2.6427606) q[2];
rz(1.8286797) q[3];
sx q[3];
rz(-2.4098318) q[3];
sx q[3];
rz(-1.7074728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6071534) q[0];
sx q[0];
rz(-0.3648912) q[0];
sx q[0];
rz(0.69951192) q[0];
rz(-2.6761159) q[1];
sx q[1];
rz(-0.87124467) q[1];
sx q[1];
rz(2.2043998) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58317157) q[0];
sx q[0];
rz(-1.4025549) q[0];
sx q[0];
rz(-2.9267163) q[0];
x q[1];
rz(-1.466339) q[2];
sx q[2];
rz(-1.8091996) q[2];
sx q[2];
rz(1.6676211) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73884642) q[1];
sx q[1];
rz(-0.35772309) q[1];
sx q[1];
rz(-0.92836942) q[1];
x q[2];
rz(1.2770416) q[3];
sx q[3];
rz(-1.2806727) q[3];
sx q[3];
rz(-1.4690831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.297544) q[2];
sx q[2];
rz(-1.301845) q[2];
sx q[2];
rz(-2.76827) q[2];
rz(1.1897872) q[3];
sx q[3];
rz(-2.6326284) q[3];
sx q[3];
rz(2.6313307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74476403) q[0];
sx q[0];
rz(-1.0823534) q[0];
sx q[0];
rz(-2.3936791) q[0];
rz(0.76639908) q[1];
sx q[1];
rz(-2.8728569) q[1];
sx q[1];
rz(0.0029729923) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9342614) q[0];
sx q[0];
rz(-0.89326477) q[0];
sx q[0];
rz(-0.24768655) q[0];
rz(-pi) q[1];
rz(-1.4998798) q[2];
sx q[2];
rz(-1.5902963) q[2];
sx q[2];
rz(-1.2649346) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.21196762) q[1];
sx q[1];
rz(-0.16008437) q[1];
sx q[1];
rz(1.2950241) q[1];
rz(2.3467196) q[3];
sx q[3];
rz(-0.27471457) q[3];
sx q[3];
rz(0.87546722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.11091867) q[2];
sx q[2];
rz(-1.7577533) q[2];
sx q[2];
rz(-2.0745011) q[2];
rz(-0.083960697) q[3];
sx q[3];
rz(-0.48560086) q[3];
sx q[3];
rz(0.77897227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.344051) q[0];
sx q[0];
rz(-2.2028956) q[0];
sx q[0];
rz(0.10636605) q[0];
rz(-0.97995177) q[1];
sx q[1];
rz(-1.6500902) q[1];
sx q[1];
rz(-0.76593691) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2625807) q[0];
sx q[0];
rz(-2.0011138) q[0];
sx q[0];
rz(2.5613214) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7226549) q[2];
sx q[2];
rz(-2.9235002) q[2];
sx q[2];
rz(-2.5401087) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.76141833) q[1];
sx q[1];
rz(-1.5360218) q[1];
sx q[1];
rz(2.4654287) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5745138) q[3];
sx q[3];
rz(-1.5848985) q[3];
sx q[3];
rz(-2.5013347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6672259) q[2];
sx q[2];
rz(-0.52046481) q[2];
sx q[2];
rz(2.4412947) q[2];
rz(0.66655603) q[3];
sx q[3];
rz(-1.3349814) q[3];
sx q[3];
rz(2.0535645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4661082) q[0];
sx q[0];
rz(-1.0270783) q[0];
sx q[0];
rz(-1.745537) q[0];
rz(-1.667977) q[1];
sx q[1];
rz(-1.2584078) q[1];
sx q[1];
rz(2.4748763) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4748866) q[0];
sx q[0];
rz(-2.1726554) q[0];
sx q[0];
rz(2.3573897) q[0];
rz(-pi) q[1];
x q[1];
rz(2.780203) q[2];
sx q[2];
rz(-0.69212428) q[2];
sx q[2];
rz(0.036594242) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.085070193) q[1];
sx q[1];
rz(-1.7016954) q[1];
sx q[1];
rz(-0.31899778) q[1];
x q[2];
rz(-2.5693043) q[3];
sx q[3];
rz(-1.0977355) q[3];
sx q[3];
rz(0.285404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0746158) q[2];
sx q[2];
rz(-2.4846027) q[2];
sx q[2];
rz(-0.89078772) q[2];
rz(0.43205076) q[3];
sx q[3];
rz(-1.0703577) q[3];
sx q[3];
rz(2.3945358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4074832) q[0];
sx q[0];
rz(-1.4984087) q[0];
sx q[0];
rz(1.8468504) q[0];
rz(1.8941849) q[1];
sx q[1];
rz(-2.4220962) q[1];
sx q[1];
rz(-1.9069506) q[1];
rz(-1.1449849) q[2];
sx q[2];
rz(-0.37526423) q[2];
sx q[2];
rz(3.0107977) q[2];
rz(0.4890781) q[3];
sx q[3];
rz(-1.5413956) q[3];
sx q[3];
rz(1.2914381) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
