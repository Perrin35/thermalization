OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0796354) q[0];
sx q[0];
rz(-2.893653) q[0];
sx q[0];
rz(0.82013446) q[0];
rz(0.040854383) q[1];
sx q[1];
rz(3.9305384) q[1];
sx q[1];
rz(9.5126704) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1782921) q[0];
sx q[0];
rz(-2.0110197) q[0];
sx q[0];
rz(-2.630169) q[0];
rz(-pi) q[1];
x q[1];
rz(2.232993) q[2];
sx q[2];
rz(-0.63630051) q[2];
sx q[2];
rz(-1.5719906) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.80884113) q[1];
sx q[1];
rz(-2.07507) q[1];
sx q[1];
rz(-1.5732952) q[1];
rz(-pi) q[2];
rz(1.8294931) q[3];
sx q[3];
rz(-1.5764578) q[3];
sx q[3];
rz(-2.2744846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0007881) q[2];
sx q[2];
rz(-1.452383) q[2];
sx q[2];
rz(0.56935707) q[2];
rz(1.5287483) q[3];
sx q[3];
rz(-0.6362392) q[3];
sx q[3];
rz(-1.7830085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59812087) q[0];
sx q[0];
rz(-2.9765029) q[0];
sx q[0];
rz(0.55364451) q[0];
rz(-1.2373295) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(-1.9083317) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2376018) q[0];
sx q[0];
rz(-1.0407789) q[0];
sx q[0];
rz(0.7941829) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1594752) q[2];
sx q[2];
rz(-2.4607686) q[2];
sx q[2];
rz(-2.8603539) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.065029) q[1];
sx q[1];
rz(-2.934888) q[1];
sx q[1];
rz(-1.4935342) q[1];
rz(-1.8096829) q[3];
sx q[3];
rz(-2.2242862) q[3];
sx q[3];
rz(1.4373923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1374986) q[2];
sx q[2];
rz(-1.6823021) q[2];
sx q[2];
rz(-3.1393576) q[2];
rz(0.8301174) q[3];
sx q[3];
rz(-2.3486962) q[3];
sx q[3];
rz(-1.8921651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8925791) q[0];
sx q[0];
rz(-2.5254288) q[0];
sx q[0];
rz(2.2900443) q[0];
rz(0.7710723) q[1];
sx q[1];
rz(-1.4665736) q[1];
sx q[1];
rz(0.071203701) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.035633798) q[0];
sx q[0];
rz(-1.7674315) q[0];
sx q[0];
rz(1.3282302) q[0];
x q[1];
rz(-2.0257646) q[2];
sx q[2];
rz(-1.3340545) q[2];
sx q[2];
rz(-0.67982212) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0752807) q[1];
sx q[1];
rz(-1.4367141) q[1];
sx q[1];
rz(-0.65384298) q[1];
rz(2.2749388) q[3];
sx q[3];
rz(-0.94368499) q[3];
sx q[3];
rz(1.8303527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8292134) q[2];
sx q[2];
rz(-2.2177314) q[2];
sx q[2];
rz(-0.52345792) q[2];
rz(-0.33189014) q[3];
sx q[3];
rz(-0.060398014) q[3];
sx q[3];
rz(0.50326842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7317384) q[0];
sx q[0];
rz(-0.24895746) q[0];
sx q[0];
rz(0.82558924) q[0];
rz(1.1666974) q[1];
sx q[1];
rz(-2.2749133) q[1];
sx q[1];
rz(1.694214) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4226845) q[0];
sx q[0];
rz(-1.8005162) q[0];
sx q[0];
rz(0.16750383) q[0];
rz(-pi) q[1];
rz(-1.6214192) q[2];
sx q[2];
rz(-2.9630337) q[2];
sx q[2];
rz(0.52390097) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.73562685) q[1];
sx q[1];
rz(-1.7298797) q[1];
sx q[1];
rz(2.2707978) q[1];
rz(-1.8131687) q[3];
sx q[3];
rz(-1.216785) q[3];
sx q[3];
rz(0.44204516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6772785) q[2];
sx q[2];
rz(-1.7769287) q[2];
sx q[2];
rz(-0.14492598) q[2];
rz(1.011301) q[3];
sx q[3];
rz(-0.62711182) q[3];
sx q[3];
rz(-3.1269126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.99575627) q[0];
sx q[0];
rz(-0.053481426) q[0];
sx q[0];
rz(-2.4575535) q[0];
rz(-1.9513291) q[1];
sx q[1];
rz(-0.66527706) q[1];
sx q[1];
rz(-1.9285944) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0096668) q[0];
sx q[0];
rz(-1.3103232) q[0];
sx q[0];
rz(0.38513222) q[0];
x q[1];
rz(2.2239387) q[2];
sx q[2];
rz(-2.0954164) q[2];
sx q[2];
rz(-1.4337002) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.1274384) q[1];
sx q[1];
rz(-1.6886097) q[1];
sx q[1];
rz(2.884306) q[1];
rz(-pi) q[2];
rz(-1.3092242) q[3];
sx q[3];
rz(-1.3993169) q[3];
sx q[3];
rz(1.8985626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.50903901) q[2];
sx q[2];
rz(-1.2712487) q[2];
sx q[2];
rz(0.63009134) q[2];
rz(-2.5693494) q[3];
sx q[3];
rz(-0.64544353) q[3];
sx q[3];
rz(-0.8852638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0260139) q[0];
sx q[0];
rz(-2.2981839) q[0];
sx q[0];
rz(0.28636006) q[0];
rz(-2.5419366) q[1];
sx q[1];
rz(-1.1279794) q[1];
sx q[1];
rz(2.5568331) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6385348) q[0];
sx q[0];
rz(-0.98031822) q[0];
sx q[0];
rz(-1.9496586) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61442394) q[2];
sx q[2];
rz(-1.5493869) q[2];
sx q[2];
rz(-1.1504088) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.11440052) q[1];
sx q[1];
rz(-0.90585867) q[1];
sx q[1];
rz(2.7726638) q[1];
rz(2.9752475) q[3];
sx q[3];
rz(-1.6206985) q[3];
sx q[3];
rz(-2.9970616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.16453234) q[2];
sx q[2];
rz(-0.88447276) q[2];
sx q[2];
rz(-1.6284846) q[2];
rz(-0.96308723) q[3];
sx q[3];
rz(-1.3685127) q[3];
sx q[3];
rz(-2.517038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0528089) q[0];
sx q[0];
rz(-1.4316906) q[0];
sx q[0];
rz(2.5928296) q[0];
rz(0.051963003) q[1];
sx q[1];
rz(-0.49182645) q[1];
sx q[1];
rz(0.60639492) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.081189922) q[0];
sx q[0];
rz(-2.5572544) q[0];
sx q[0];
rz(-2.0085745) q[0];
x q[1];
rz(2.0386001) q[2];
sx q[2];
rz(-2.4228547) q[2];
sx q[2];
rz(-1.7407835) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.28601521) q[1];
sx q[1];
rz(-0.77076036) q[1];
sx q[1];
rz(-1.4090178) q[1];
x q[2];
rz(-1.9803195) q[3];
sx q[3];
rz(-1.2425353) q[3];
sx q[3];
rz(1.3271774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9397883) q[2];
sx q[2];
rz(-2.1494892) q[2];
sx q[2];
rz(-2.6331804) q[2];
rz(-1.1172179) q[3];
sx q[3];
rz(-2.9682187) q[3];
sx q[3];
rz(1.6569051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8023119) q[0];
sx q[0];
rz(-2.7957714) q[0];
sx q[0];
rz(-0.75396496) q[0];
rz(-0.94999653) q[1];
sx q[1];
rz(-1.4818622) q[1];
sx q[1];
rz(2.9139013) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65093016) q[0];
sx q[0];
rz(-2.7224446) q[0];
sx q[0];
rz(1.3952257) q[0];
x q[1];
rz(-2.6954306) q[2];
sx q[2];
rz(-1.9036306) q[2];
sx q[2];
rz(1.7762426) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7714872) q[1];
sx q[1];
rz(-2.5574554) q[1];
sx q[1];
rz(2.7258956) q[1];
x q[2];
rz(-3.0759096) q[3];
sx q[3];
rz(-2.5816397) q[3];
sx q[3];
rz(-0.63510676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.078538744) q[2];
sx q[2];
rz(-2.9324014) q[2];
sx q[2];
rz(0.83855808) q[2];
rz(2.7770384) q[3];
sx q[3];
rz(-1.3804599) q[3];
sx q[3];
rz(2.3013039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.7446328) q[0];
sx q[0];
rz(-1.7663706) q[0];
sx q[0];
rz(-0.80378419) q[0];
rz(-1.0546168) q[1];
sx q[1];
rz(-0.63916731) q[1];
sx q[1];
rz(-1.1484336) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9371985) q[0];
sx q[0];
rz(-3.0657401) q[0];
sx q[0];
rz(1.2501459) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0954275) q[2];
sx q[2];
rz(-1.2534007) q[2];
sx q[2];
rz(1.5517) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.25363898) q[1];
sx q[1];
rz(-1.1432853) q[1];
sx q[1];
rz(-1.0484139) q[1];
rz(-pi) q[2];
x q[2];
rz(0.095454772) q[3];
sx q[3];
rz(-2.7881552) q[3];
sx q[3];
rz(0.27775882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.049872963) q[2];
sx q[2];
rz(-1.0174948) q[2];
sx q[2];
rz(2.3633374) q[2];
rz(2.9459279) q[3];
sx q[3];
rz(-1.0396495) q[3];
sx q[3];
rz(-1.0218609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2739928) q[0];
sx q[0];
rz(-0.99853981) q[0];
sx q[0];
rz(-2.8883949) q[0];
rz(-2.4018438) q[1];
sx q[1];
rz(-2.4052129) q[1];
sx q[1];
rz(2.443312) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88903504) q[0];
sx q[0];
rz(-1.4989984) q[0];
sx q[0];
rz(0.30307146) q[0];
x q[1];
rz(3.1168078) q[2];
sx q[2];
rz(-1.6449271) q[2];
sx q[2];
rz(-2.5075846) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.56863927) q[1];
sx q[1];
rz(-1.1146953) q[1];
sx q[1];
rz(-0.24104636) q[1];
x q[2];
rz(-1.0032759) q[3];
sx q[3];
rz(-1.365981) q[3];
sx q[3];
rz(1.8703465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11761052) q[2];
sx q[2];
rz(-0.88279804) q[2];
sx q[2];
rz(2.3013766) q[2];
rz(2.501287) q[3];
sx q[3];
rz(-0.87602031) q[3];
sx q[3];
rz(1.1673814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2587851) q[0];
sx q[0];
rz(-0.85492815) q[0];
sx q[0];
rz(1.7229765) q[0];
rz(-3.1124658) q[1];
sx q[1];
rz(-0.11321414) q[1];
sx q[1];
rz(-1.3197457) q[1];
rz(-0.48568934) q[2];
sx q[2];
rz(-2.4808241) q[2];
sx q[2];
rz(-0.11447699) q[2];
rz(-1.7520262) q[3];
sx q[3];
rz(-2.443234) q[3];
sx q[3];
rz(0.20886226) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
