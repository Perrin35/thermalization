OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.68706694) q[0];
sx q[0];
rz(-1.878976) q[0];
sx q[0];
rz(-0.73448056) q[0];
rz(1.6425411) q[1];
sx q[1];
rz(1.8416815) q[1];
sx q[1];
rz(11.587974) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9578732) q[0];
sx q[0];
rz(-1.5786202) q[0];
sx q[0];
rz(1.70723) q[0];
x q[1];
rz(-2.5402563) q[2];
sx q[2];
rz(-2.3365855) q[2];
sx q[2];
rz(2.7108101) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9339188) q[1];
sx q[1];
rz(-1.7130995) q[1];
sx q[1];
rz(1.2486904) q[1];
x q[2];
rz(0.079732883) q[3];
sx q[3];
rz(-2.0679722) q[3];
sx q[3];
rz(-0.91547188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.80028278) q[2];
sx q[2];
rz(-1.6494696) q[2];
sx q[2];
rz(2.4260803) q[2];
rz(1.7320775) q[3];
sx q[3];
rz(-2.2655497) q[3];
sx q[3];
rz(1.5040262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3187123) q[0];
sx q[0];
rz(-2.5647698) q[0];
sx q[0];
rz(0.034828287) q[0];
rz(-2.4233129) q[1];
sx q[1];
rz(-1.4052582) q[1];
sx q[1];
rz(1.1911596) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1930192) q[0];
sx q[0];
rz(-2.5002242) q[0];
sx q[0];
rz(-1.8383584) q[0];
x q[1];
rz(1.7188363) q[2];
sx q[2];
rz(-2.2477373) q[2];
sx q[2];
rz(-1.9612775) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4264448) q[1];
sx q[1];
rz(-2.1472125) q[1];
sx q[1];
rz(1.2395241) q[1];
rz(-pi) q[2];
rz(2.0091811) q[3];
sx q[3];
rz(-1.7436641) q[3];
sx q[3];
rz(2.5684484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8027773) q[2];
sx q[2];
rz(-1.85314) q[2];
sx q[2];
rz(0.046099376) q[2];
rz(2.3378546) q[3];
sx q[3];
rz(-1.2396783) q[3];
sx q[3];
rz(2.5900335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7130724) q[0];
sx q[0];
rz(-1.5794733) q[0];
sx q[0];
rz(2.5285316) q[0];
rz(0.26519457) q[1];
sx q[1];
rz(-1.3399597) q[1];
sx q[1];
rz(-3.0934966) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10079028) q[0];
sx q[0];
rz(-1.5682966) q[0];
sx q[0];
rz(-2.0526171) q[0];
x q[1];
rz(-2.0916284) q[2];
sx q[2];
rz(-2.0893059) q[2];
sx q[2];
rz(1.1074378) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4029419) q[1];
sx q[1];
rz(-1.5183107) q[1];
sx q[1];
rz(1.5110635) q[1];
rz(-pi) q[2];
rz(-0.35086029) q[3];
sx q[3];
rz(-1.8985629) q[3];
sx q[3];
rz(1.8192489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.6951822) q[2];
sx q[2];
rz(-1.1626652) q[2];
sx q[2];
rz(-2.7282558) q[2];
rz(-0.6684331) q[3];
sx q[3];
rz(-1.3276525) q[3];
sx q[3];
rz(-1.7641164) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9628825) q[0];
sx q[0];
rz(-0.2186192) q[0];
sx q[0];
rz(2.9982153) q[0];
rz(-0.42477056) q[1];
sx q[1];
rz(-1.2062585) q[1];
sx q[1];
rz(0.43407789) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052267678) q[0];
sx q[0];
rz(-0.10073034) q[0];
sx q[0];
rz(1.9585376) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9888568) q[2];
sx q[2];
rz(-0.65780168) q[2];
sx q[2];
rz(1.3721823) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39555711) q[1];
sx q[1];
rz(-1.0537123) q[1];
sx q[1];
rz(-0.16090572) q[1];
x q[2];
rz(1.7884718) q[3];
sx q[3];
rz(-0.77131144) q[3];
sx q[3];
rz(1.9096489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8726337) q[2];
sx q[2];
rz(-2.3838398) q[2];
sx q[2];
rz(2.7244275) q[2];
rz(-2.4560691) q[3];
sx q[3];
rz(-1.7020107) q[3];
sx q[3];
rz(-0.050366966) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79614574) q[0];
sx q[0];
rz(-2.8247927) q[0];
sx q[0];
rz(1.5078804) q[0];
rz(-0.66868526) q[1];
sx q[1];
rz(-1.94328) q[1];
sx q[1];
rz(-2.0100458) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39172669) q[0];
sx q[0];
rz(-1.5380812) q[0];
sx q[0];
rz(0.97264887) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8077233) q[2];
sx q[2];
rz(-1.0354932) q[2];
sx q[2];
rz(-1.6795422) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.9100489) q[1];
sx q[1];
rz(-2.0374415) q[1];
sx q[1];
rz(-1.779516) q[1];
rz(0.11891811) q[3];
sx q[3];
rz(-1.44373) q[3];
sx q[3];
rz(2.6112664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9097627) q[2];
sx q[2];
rz(-2.6015687) q[2];
sx q[2];
rz(0.46169272) q[2];
rz(0.26990226) q[3];
sx q[3];
rz(-1.2609127) q[3];
sx q[3];
rz(-0.65129483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5676024) q[0];
sx q[0];
rz(-0.42581588) q[0];
sx q[0];
rz(2.0800166) q[0];
rz(2.4827982) q[1];
sx q[1];
rz(-1.7196722) q[1];
sx q[1];
rz(-2.7366791) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18428567) q[0];
sx q[0];
rz(-0.71249639) q[0];
sx q[0];
rz(-0.16100968) q[0];
rz(2.5091841) q[2];
sx q[2];
rz(-2.5002865) q[2];
sx q[2];
rz(-0.270688) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9820653) q[1];
sx q[1];
rz(-1.7100466) q[1];
sx q[1];
rz(2.0811446) q[1];
rz(1.9642716) q[3];
sx q[3];
rz(-0.67613908) q[3];
sx q[3];
rz(-2.2666988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.095526516) q[2];
sx q[2];
rz(-2.1370856) q[2];
sx q[2];
rz(1.7556165) q[2];
rz(0.69685495) q[3];
sx q[3];
rz(-0.98074061) q[3];
sx q[3];
rz(-2.3314886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.54414576) q[0];
sx q[0];
rz(-0.61328855) q[0];
sx q[0];
rz(0.59462732) q[0];
rz(2.0095297) q[1];
sx q[1];
rz(-1.4040399) q[1];
sx q[1];
rz(-0.51116699) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5385732) q[0];
sx q[0];
rz(-1.5165189) q[0];
sx q[0];
rz(0.57247573) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2726502) q[2];
sx q[2];
rz(-0.41837439) q[2];
sx q[2];
rz(-1.6623896) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0446544) q[1];
sx q[1];
rz(-1.6676098) q[1];
sx q[1];
rz(-1.8332048) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4092402) q[3];
sx q[3];
rz(-1.9480289) q[3];
sx q[3];
rz(-2.4767552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1749997) q[2];
sx q[2];
rz(-0.24093691) q[2];
sx q[2];
rz(0.28755507) q[2];
rz(1.1047085) q[3];
sx q[3];
rz(-1.46773) q[3];
sx q[3];
rz(-0.12565676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48968807) q[0];
sx q[0];
rz(-2.7629485) q[0];
sx q[0];
rz(-2.4878159) q[0];
rz(-2.6405624) q[1];
sx q[1];
rz(-2.4153695) q[1];
sx q[1];
rz(2.3562866) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57779271) q[0];
sx q[0];
rz(-2.1350088) q[0];
sx q[0];
rz(1.7336506) q[0];
x q[1];
rz(-2.4926033) q[2];
sx q[2];
rz(-2.4675998) q[2];
sx q[2];
rz(1.4907995) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0846906) q[1];
sx q[1];
rz(-0.71699079) q[1];
sx q[1];
rz(-0.3221953) q[1];
rz(-pi) q[2];
rz(-2.9124746) q[3];
sx q[3];
rz(-1.4909298) q[3];
sx q[3];
rz(-0.95019482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.952897) q[2];
sx q[2];
rz(-1.5263824) q[2];
sx q[2];
rz(-0.3696028) q[2];
rz(2.8772307) q[3];
sx q[3];
rz(-1.8601469) q[3];
sx q[3];
rz(3.1407691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13067506) q[0];
sx q[0];
rz(-2.4322746) q[0];
sx q[0];
rz(1.5189019) q[0];
rz(1.9021289) q[1];
sx q[1];
rz(-2.0294956) q[1];
sx q[1];
rz(2.1942203) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1993123) q[0];
sx q[0];
rz(-0.40422842) q[0];
sx q[0];
rz(1.4154427) q[0];
x q[1];
rz(2.9765509) q[2];
sx q[2];
rz(-0.36569917) q[2];
sx q[2];
rz(-0.62830454) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.92129642) q[1];
sx q[1];
rz(-1.4054757) q[1];
sx q[1];
rz(-2.9784747) q[1];
rz(-pi) q[2];
rz(1.4072284) q[3];
sx q[3];
rz(-1.6298098) q[3];
sx q[3];
rz(1.407377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.47306791) q[2];
sx q[2];
rz(-0.99978414) q[2];
sx q[2];
rz(-1.1119941) q[2];
rz(-0.2019349) q[3];
sx q[3];
rz(-0.58378059) q[3];
sx q[3];
rz(0.37566617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9428381) q[0];
sx q[0];
rz(-0.18809479) q[0];
sx q[0];
rz(1.4659708) q[0];
rz(-1.6000043) q[1];
sx q[1];
rz(-1.8406638) q[1];
sx q[1];
rz(0.25513729) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32624558) q[0];
sx q[0];
rz(-1.4370942) q[0];
sx q[0];
rz(-3.1241425) q[0];
x q[1];
rz(-0.77346797) q[2];
sx q[2];
rz(-1.3313401) q[2];
sx q[2];
rz(2.4635893) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0855116) q[1];
sx q[1];
rz(-1.8485118) q[1];
sx q[1];
rz(-3.059832) q[1];
rz(-pi) q[2];
rz(0.76948036) q[3];
sx q[3];
rz(-1.7090685) q[3];
sx q[3];
rz(-1.4797803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.61350358) q[2];
sx q[2];
rz(-1.1897503) q[2];
sx q[2];
rz(2.7733754) q[2];
rz(0.28299371) q[3];
sx q[3];
rz(-0.26809433) q[3];
sx q[3];
rz(0.32576573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014864347) q[0];
sx q[0];
rz(-2.2751502) q[0];
sx q[0];
rz(-1.0979102) q[0];
rz(2.6425843) q[1];
sx q[1];
rz(-1.9022763) q[1];
sx q[1];
rz(0.3872445) q[1];
rz(2.8088612) q[2];
sx q[2];
rz(-1.6920148) q[2];
sx q[2];
rz(0.68670338) q[2];
rz(-3.079703) q[3];
sx q[3];
rz(-1.7714942) q[3];
sx q[3];
rz(-1.9006806) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
