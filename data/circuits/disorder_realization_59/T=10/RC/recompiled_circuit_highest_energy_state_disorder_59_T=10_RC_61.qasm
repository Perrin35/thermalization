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
rz(1.264313) q[0];
sx q[0];
rz(-2.8180583) q[0];
sx q[0];
rz(-2.6927595) q[0];
rz(-1.1558865) q[1];
sx q[1];
rz(-1.7855676) q[1];
sx q[1];
rz(1.9670271) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1313707) q[0];
sx q[0];
rz(-1.0099942) q[0];
sx q[0];
rz(1.9813374) q[0];
x q[1];
rz(-2.817906) q[2];
sx q[2];
rz(-1.5059221) q[2];
sx q[2];
rz(1.5144358) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5331363) q[1];
sx q[1];
rz(-1.605833) q[1];
sx q[1];
rz(2.1916981) q[1];
x q[2];
rz(1.6714736) q[3];
sx q[3];
rz(-1.0591456) q[3];
sx q[3];
rz(2.4477521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55149469) q[2];
sx q[2];
rz(-1.3000725) q[2];
sx q[2];
rz(-0.91280118) q[2];
rz(-1.8707976) q[3];
sx q[3];
rz(-2.1015034) q[3];
sx q[3];
rz(2.2916268) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1368644) q[0];
sx q[0];
rz(-0.57924634) q[0];
sx q[0];
rz(0.42214033) q[0];
rz(-0.36960754) q[1];
sx q[1];
rz(-2.044675) q[1];
sx q[1];
rz(-0.94013989) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63623171) q[0];
sx q[0];
rz(-2.8121037) q[0];
sx q[0];
rz(-1.5630653) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1040078) q[2];
sx q[2];
rz(-1.4681446) q[2];
sx q[2];
rz(0.29666049) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.2804037) q[1];
sx q[1];
rz(-2.197816) q[1];
sx q[1];
rz(-0.95071937) q[1];
rz(-pi) q[2];
rz(0.89408447) q[3];
sx q[3];
rz(-1.3562804) q[3];
sx q[3];
rz(1.0971951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9422354) q[2];
sx q[2];
rz(-0.67023674) q[2];
sx q[2];
rz(-2.2352236) q[2];
rz(0.97936112) q[3];
sx q[3];
rz(-2.7370079) q[3];
sx q[3];
rz(1.8766859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1449428) q[0];
sx q[0];
rz(-3.0746089) q[0];
sx q[0];
rz(1.2445194) q[0];
rz(0.38885802) q[1];
sx q[1];
rz(-1.2797979) q[1];
sx q[1];
rz(-0.32274524) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6038216) q[0];
sx q[0];
rz(-2.4706998) q[0];
sx q[0];
rz(0.18108271) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1305442) q[2];
sx q[2];
rz(-1.3072392) q[2];
sx q[2];
rz(-1.6853058) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6420329) q[1];
sx q[1];
rz(-0.82553464) q[1];
sx q[1];
rz(-0.31742974) q[1];
rz(-1.6204319) q[3];
sx q[3];
rz(-1.2200886) q[3];
sx q[3];
rz(-1.3405104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2286223) q[2];
sx q[2];
rz(-1.7614438) q[2];
sx q[2];
rz(1.3529533) q[2];
rz(-2.8690423) q[3];
sx q[3];
rz(-1.8504146) q[3];
sx q[3];
rz(-1.6779617) q[3];
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
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54101855) q[0];
sx q[0];
rz(-0.32222846) q[0];
sx q[0];
rz(-1.2433276) q[0];
rz(2.2728032) q[1];
sx q[1];
rz(-0.96005762) q[1];
sx q[1];
rz(2.0713846) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7420827) q[0];
sx q[0];
rz(-2.1993981) q[0];
sx q[0];
rz(2.5058772) q[0];
rz(-2.0212428) q[2];
sx q[2];
rz(-2.3131158) q[2];
sx q[2];
rz(2.7591005) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3096042) q[1];
sx q[1];
rz(-1.8924531) q[1];
sx q[1];
rz(-3.020806) q[1];
rz(-pi) q[2];
rz(0.34262212) q[3];
sx q[3];
rz(-1.3903769) q[3];
sx q[3];
rz(2.8460549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.59226817) q[2];
sx q[2];
rz(-2.3531239) q[2];
sx q[2];
rz(-2.2198524) q[2];
rz(-2.8978469) q[3];
sx q[3];
rz(-1.4344401) q[3];
sx q[3];
rz(-0.29269472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9859966) q[0];
sx q[0];
rz(-1.3256185) q[0];
sx q[0];
rz(2.6407114) q[0];
rz(-2.0241375) q[1];
sx q[1];
rz(-1.8084348) q[1];
sx q[1];
rz(-2.8270328) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4582493) q[0];
sx q[0];
rz(-1.9106081) q[0];
sx q[0];
rz(0.37796867) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0044745046) q[2];
sx q[2];
rz(-1.5273849) q[2];
sx q[2];
rz(2.2501066) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6198338) q[1];
sx q[1];
rz(-1.17621) q[1];
sx q[1];
rz(2.16741) q[1];
rz(-pi) q[2];
rz(-0.1184095) q[3];
sx q[3];
rz(-1.1573912) q[3];
sx q[3];
rz(1.4865231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.306281) q[2];
sx q[2];
rz(-2.0302782) q[2];
sx q[2];
rz(1.296898) q[2];
rz(-0.48526192) q[3];
sx q[3];
rz(-1.5596215) q[3];
sx q[3];
rz(-2.0280793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49419633) q[0];
sx q[0];
rz(-1.2911456) q[0];
sx q[0];
rz(-3.0418292) q[0];
rz(-1.8708723) q[1];
sx q[1];
rz(-2.2427509) q[1];
sx q[1];
rz(1.3894003) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0770806) q[0];
sx q[0];
rz(-1.837366) q[0];
sx q[0];
rz(1.7414016) q[0];
rz(-0.17398457) q[2];
sx q[2];
rz(-1.4672888) q[2];
sx q[2];
rz(0.058920842) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.8627117) q[1];
sx q[1];
rz(-0.91966719) q[1];
sx q[1];
rz(2.34823) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12677315) q[3];
sx q[3];
rz(-0.78713464) q[3];
sx q[3];
rz(-1.2490369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6483267) q[2];
sx q[2];
rz(-2.6706225) q[2];
sx q[2];
rz(-2.2512482) q[2];
rz(1.1912311) q[3];
sx q[3];
rz(-1.3072562) q[3];
sx q[3];
rz(0.90203917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4297727) q[0];
sx q[0];
rz(-3.0968102) q[0];
sx q[0];
rz(-0.032715948) q[0];
rz(2.6383767) q[1];
sx q[1];
rz(-1.0992173) q[1];
sx q[1];
rz(3.018766) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5779992) q[0];
sx q[0];
rz(-1.8933354) q[0];
sx q[0];
rz(-2.3911227) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5678504) q[2];
sx q[2];
rz(-2.5649568) q[2];
sx q[2];
rz(1.6078311) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.87007755) q[1];
sx q[1];
rz(-2.6598499) q[1];
sx q[1];
rz(-1.7863356) q[1];
rz(0.35251725) q[3];
sx q[3];
rz(-1.4488701) q[3];
sx q[3];
rz(-2.8059949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7922625) q[2];
sx q[2];
rz(-0.96668875) q[2];
sx q[2];
rz(1.1518504) q[2];
rz(0.88519111) q[3];
sx q[3];
rz(-1.7722292) q[3];
sx q[3];
rz(1.4256029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18871466) q[0];
sx q[0];
rz(-0.9335683) q[0];
sx q[0];
rz(2.2123912) q[0];
rz(-2.5125391) q[1];
sx q[1];
rz(-1.7864497) q[1];
sx q[1];
rz(0.74310511) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6154872) q[0];
sx q[0];
rz(-2.0172523) q[0];
sx q[0];
rz(-0.087421405) q[0];
x q[1];
rz(1.475072) q[2];
sx q[2];
rz(-2.0572955) q[2];
sx q[2];
rz(-2.7216788) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9070702) q[1];
sx q[1];
rz(-1.6369695) q[1];
sx q[1];
rz(0.41445865) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2065714) q[3];
sx q[3];
rz(-1.1639769) q[3];
sx q[3];
rz(1.1233028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5653845) q[2];
sx q[2];
rz(-2.3614063) q[2];
sx q[2];
rz(0.72594491) q[2];
rz(-2.7706326) q[3];
sx q[3];
rz(-1.5981263) q[3];
sx q[3];
rz(0.36271873) q[3];
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
rz(-2.2711656) q[0];
sx q[0];
rz(-0.83036625) q[0];
sx q[0];
rz(-0.56394947) q[0];
rz(-1.14934) q[1];
sx q[1];
rz(-2.1795858) q[1];
sx q[1];
rz(-2.3990778) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71096984) q[0];
sx q[0];
rz(-1.1046243) q[0];
sx q[0];
rz(-2.6197893) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41461035) q[2];
sx q[2];
rz(-1.3002204) q[2];
sx q[2];
rz(2.0913578) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0997252) q[1];
sx q[1];
rz(-2.0085287) q[1];
sx q[1];
rz(-1.5220286) q[1];
x q[2];
rz(2.5768859) q[3];
sx q[3];
rz(-0.93834025) q[3];
sx q[3];
rz(-2.6868827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.62770647) q[2];
sx q[2];
rz(-2.0676925) q[2];
sx q[2];
rz(0.34454301) q[2];
rz(-0.44376093) q[3];
sx q[3];
rz(-1.0756451) q[3];
sx q[3];
rz(-1.3226604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3798856) q[0];
sx q[0];
rz(-2.8590617) q[0];
sx q[0];
rz(-0.66201061) q[0];
rz(-2.8061891) q[1];
sx q[1];
rz(-1.4619724) q[1];
sx q[1];
rz(-0.9368771) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4136018) q[0];
sx q[0];
rz(-1.2112036) q[0];
sx q[0];
rz(0.92443621) q[0];
rz(-1.4858021) q[2];
sx q[2];
rz(-2.8704145) q[2];
sx q[2];
rz(2.8176339) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.42439991) q[1];
sx q[1];
rz(-2.6584627) q[1];
sx q[1];
rz(0.65349726) q[1];
x q[2];
rz(-2.9659445) q[3];
sx q[3];
rz(-1.3743168) q[3];
sx q[3];
rz(0.59284569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.718049) q[2];
sx q[2];
rz(-0.23423883) q[2];
sx q[2];
rz(2.6424109) q[2];
rz(-1.6792123) q[3];
sx q[3];
rz(-1.1464109) q[3];
sx q[3];
rz(-1.6818989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.725631) q[0];
sx q[0];
rz(-1.220663) q[0];
sx q[0];
rz(2.2882373) q[0];
rz(2.536643) q[1];
sx q[1];
rz(-1.1440944) q[1];
sx q[1];
rz(0.49857421) q[1];
rz(-1.4251322) q[2];
sx q[2];
rz(-2.1389516) q[2];
sx q[2];
rz(0.4898796) q[2];
rz(-0.50157401) q[3];
sx q[3];
rz(-1.4091103) q[3];
sx q[3];
rz(1.4179358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
