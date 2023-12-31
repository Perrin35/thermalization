OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2621736) q[0];
sx q[0];
rz(-1.7466495) q[0];
sx q[0];
rz(3.1403132) q[0];
rz(1.4446422) q[1];
sx q[1];
rz(-1.1029707) q[1];
sx q[1];
rz(2.3666518) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89844184) q[0];
sx q[0];
rz(-1.4282707) q[0];
sx q[0];
rz(-1.5205163) q[0];
rz(-pi) q[1];
rz(-0.84471976) q[2];
sx q[2];
rz(-1.2436927) q[2];
sx q[2];
rz(0.5288045) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.70667627) q[1];
sx q[1];
rz(-1.7129363) q[1];
sx q[1];
rz(0.040277004) q[1];
x q[2];
rz(0.48304708) q[3];
sx q[3];
rz(-2.8222198) q[3];
sx q[3];
rz(-1.7628302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0960192) q[2];
sx q[2];
rz(-2.0298268) q[2];
sx q[2];
rz(-1.1958896) q[2];
rz(1.9879509) q[3];
sx q[3];
rz(-2.3524645) q[3];
sx q[3];
rz(-1.7529863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4202704) q[0];
sx q[0];
rz(-2.7852311) q[0];
sx q[0];
rz(1.2715682) q[0];
rz(-2.0416073) q[1];
sx q[1];
rz(-1.1106691) q[1];
sx q[1];
rz(-1.3756479) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028582024) q[0];
sx q[0];
rz(-2.7227289) q[0];
sx q[0];
rz(-0.21582614) q[0];
rz(3.0963418) q[2];
sx q[2];
rz(-1.3327193) q[2];
sx q[2];
rz(2.0766052) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8010506) q[1];
sx q[1];
rz(-1.2060296) q[1];
sx q[1];
rz(-1.7629452) q[1];
rz(-1.0248915) q[3];
sx q[3];
rz(-2.5149269) q[3];
sx q[3];
rz(-2.6729667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9282844) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(0.48970547) q[2];
rz(-2.0080163) q[3];
sx q[3];
rz(-2.9462892) q[3];
sx q[3];
rz(1.0666696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60004822) q[0];
sx q[0];
rz(-1.2788037) q[0];
sx q[0];
rz(0.24060732) q[0];
rz(-2.799017) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(1.906357) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9440207) q[0];
sx q[0];
rz(-0.70957843) q[0];
sx q[0];
rz(-0.59528415) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58358242) q[2];
sx q[2];
rz(-2.7036813) q[2];
sx q[2];
rz(-0.32354087) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2521378) q[1];
sx q[1];
rz(-2.7733907) q[1];
sx q[1];
rz(-2.3705269) q[1];
x q[2];
rz(-0.62927411) q[3];
sx q[3];
rz(-2.4089703) q[3];
sx q[3];
rz(-2.7936558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.76413313) q[2];
sx q[2];
rz(-1.5185792) q[2];
sx q[2];
rz(1.4366478) q[2];
rz(-1.7403587) q[3];
sx q[3];
rz(-1.2652206) q[3];
sx q[3];
rz(-0.60825545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075994611) q[0];
sx q[0];
rz(-1.6125212) q[0];
sx q[0];
rz(0.80379379) q[0];
rz(-2.1919788) q[1];
sx q[1];
rz(-1.6751553) q[1];
sx q[1];
rz(0.11985699) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45194295) q[0];
sx q[0];
rz(-1.6065238) q[0];
sx q[0];
rz(0.96468754) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36977936) q[2];
sx q[2];
rz(-0.38447194) q[2];
sx q[2];
rz(2.0535) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2377492) q[1];
sx q[1];
rz(-1.3386054) q[1];
sx q[1];
rz(-0.018149908) q[1];
rz(-pi) q[2];
rz(-0.5991163) q[3];
sx q[3];
rz(-1.1154419) q[3];
sx q[3];
rz(-1.276254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5084761) q[2];
sx q[2];
rz(-1.2458331) q[2];
sx q[2];
rz(3.0692549) q[2];
rz(0.37483254) q[3];
sx q[3];
rz(-2.4852677) q[3];
sx q[3];
rz(-1.9434631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4478093) q[0];
sx q[0];
rz(-2.5700975) q[0];
sx q[0];
rz(-2.6547292) q[0];
rz(-2.411719) q[1];
sx q[1];
rz(-0.90881538) q[1];
sx q[1];
rz(-1.240085) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70532521) q[0];
sx q[0];
rz(-1.8099394) q[0];
sx q[0];
rz(-1.5572085) q[0];
rz(-pi) q[1];
rz(0.24057062) q[2];
sx q[2];
rz(-1.8261989) q[2];
sx q[2];
rz(2.0672928) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1357437) q[1];
sx q[1];
rz(-1.504335) q[1];
sx q[1];
rz(-0.60484109) q[1];
rz(-pi) q[2];
rz(0.48360444) q[3];
sx q[3];
rz(-1.0217474) q[3];
sx q[3];
rz(-2.1476114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1029677) q[2];
sx q[2];
rz(-2.2323148) q[2];
sx q[2];
rz(2.4385578) q[2];
rz(-1.4098343) q[3];
sx q[3];
rz(-1.792428) q[3];
sx q[3];
rz(2.9838802) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1773961) q[0];
sx q[0];
rz(-1.2921768) q[0];
sx q[0];
rz(-0.99037209) q[0];
rz(0.05274996) q[1];
sx q[1];
rz(-0.92664781) q[1];
sx q[1];
rz(1.4809158) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71125644) q[0];
sx q[0];
rz(-2.7529786) q[0];
sx q[0];
rz(1.7646164) q[0];
rz(-pi) q[1];
rz(-0.00077498282) q[2];
sx q[2];
rz(-1.1256071) q[2];
sx q[2];
rz(2.0489401) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.90658014) q[1];
sx q[1];
rz(-1.3624411) q[1];
sx q[1];
rz(1.033003) q[1];
rz(-pi) q[2];
rz(-1.6310286) q[3];
sx q[3];
rz(-0.67959736) q[3];
sx q[3];
rz(2.4356902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0084373077) q[2];
sx q[2];
rz(-1.6494273) q[2];
sx q[2];
rz(2.5793502) q[2];
rz(2.0810614) q[3];
sx q[3];
rz(-2.3980467) q[3];
sx q[3];
rz(0.26091584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19787191) q[0];
sx q[0];
rz(-1.3523538) q[0];
sx q[0];
rz(-1.6725756) q[0];
rz(2.127227) q[1];
sx q[1];
rz(-1.0275774) q[1];
sx q[1];
rz(-1.8168824) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4949188) q[0];
sx q[0];
rz(-3.0122628) q[0];
sx q[0];
rz(-2.9652251) q[0];
rz(-pi) q[1];
rz(2.1340738) q[2];
sx q[2];
rz(-1.129732) q[2];
sx q[2];
rz(-1.2457459) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3130256) q[1];
sx q[1];
rz(-1.9777858) q[1];
sx q[1];
rz(-1.160497) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7571194) q[3];
sx q[3];
rz(-2.1914346) q[3];
sx q[3];
rz(-1.4409161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5801195) q[2];
sx q[2];
rz(-1.7529528) q[2];
sx q[2];
rz(0.44357792) q[2];
rz(-2.1929072) q[3];
sx q[3];
rz(-1.9493999) q[3];
sx q[3];
rz(0.77478066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.7396486) q[0];
sx q[0];
rz(-2.1304603) q[0];
sx q[0];
rz(-2.6810714) q[0];
rz(0.1000239) q[1];
sx q[1];
rz(-0.99383751) q[1];
sx q[1];
rz(-1.8519648) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3562718) q[0];
sx q[0];
rz(-2.3873781) q[0];
sx q[0];
rz(-0.69099364) q[0];
x q[1];
rz(-2.3233534) q[2];
sx q[2];
rz(-2.1365676) q[2];
sx q[2];
rz(0.73275369) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3370812) q[1];
sx q[1];
rz(-1.6816499) q[1];
sx q[1];
rz(1.3827419) q[1];
rz(1.7361705) q[3];
sx q[3];
rz(-2.227042) q[3];
sx q[3];
rz(1.0429494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.148968) q[2];
sx q[2];
rz(-0.31105369) q[2];
sx q[2];
rz(-1.0827433) q[2];
rz(-3.0454214) q[3];
sx q[3];
rz(-1.7008737) q[3];
sx q[3];
rz(-1.2307897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7643395) q[0];
sx q[0];
rz(-2.9057673) q[0];
sx q[0];
rz(0.33426958) q[0];
rz(-1.9175247) q[1];
sx q[1];
rz(-1.5806438) q[1];
sx q[1];
rz(-2.8589378) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5196913) q[0];
sx q[0];
rz(-1.1150517) q[0];
sx q[0];
rz(0.68463188) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7868144) q[2];
sx q[2];
rz(-1.2470494) q[2];
sx q[2];
rz(0.87460364) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5987451) q[1];
sx q[1];
rz(-1.3951021) q[1];
sx q[1];
rz(-2.7302242) q[1];
x q[2];
rz(-1.71404) q[3];
sx q[3];
rz(-0.69056615) q[3];
sx q[3];
rz(-2.7018202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9294372) q[2];
sx q[2];
rz(-1.9781457) q[2];
sx q[2];
rz(2.4003417) q[2];
rz(-0.50179982) q[3];
sx q[3];
rz(-1.3391677) q[3];
sx q[3];
rz(0.58399502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0989477) q[0];
sx q[0];
rz(-2.2462923) q[0];
sx q[0];
rz(2.6328971) q[0];
rz(3.0264061) q[1];
sx q[1];
rz(-0.70786628) q[1];
sx q[1];
rz(2.4597816) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0081351) q[0];
sx q[0];
rz(-1.7383766) q[0];
sx q[0];
rz(2.125678) q[0];
rz(-2.290906) q[2];
sx q[2];
rz(-0.88973532) q[2];
sx q[2];
rz(-2.1249287) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0360003) q[1];
sx q[1];
rz(-0.75165527) q[1];
sx q[1];
rz(-2.4113301) q[1];
rz(-pi) q[2];
rz(-2.1308594) q[3];
sx q[3];
rz(-2.0252953) q[3];
sx q[3];
rz(-2.0927932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8455785) q[2];
sx q[2];
rz(-1.6342376) q[2];
sx q[2];
rz(-0.34109035) q[2];
rz(1.0836481) q[3];
sx q[3];
rz(-0.79569474) q[3];
sx q[3];
rz(2.9705689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1682128) q[0];
sx q[0];
rz(-1.4773049) q[0];
sx q[0];
rz(2.1422577) q[0];
rz(0.60733168) q[1];
sx q[1];
rz(-0.9098396) q[1];
sx q[1];
rz(0.20914016) q[1];
rz(-0.66424673) q[2];
sx q[2];
rz(-1.9953809) q[2];
sx q[2];
rz(-2.8473163) q[2];
rz(0.38469436) q[3];
sx q[3];
rz(-1.3833429) q[3];
sx q[3];
rz(-0.27237567) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
