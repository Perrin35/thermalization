OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.126625210046768) q[0];
sx q[0];
rz(1.57464829285676) q[0];
sx q[0];
rz(9.94045416115924) q[0];
rz(0.2281783670187) q[1];
sx q[1];
rz(5.53622666199739) q[1];
sx q[1];
rz(9.85103917717143) q[1];
cx q[1],q[0];
rz(0.928620636463165) q[0];
sx q[0];
rz(5.47336211998994) q[0];
sx q[0];
rz(8.09311399459049) q[0];
rz(1.078329205513) q[2];
sx q[2];
rz(3.70028457243974) q[2];
sx q[2];
rz(9.87148059009715) q[2];
cx q[2],q[1];
rz(-1.96102404594421) q[1];
sx q[1];
rz(3.74605444272096) q[1];
sx q[1];
rz(11.0971259832303) q[1];
rz(-2.55741381645203) q[3];
sx q[3];
rz(4.98625031312043) q[3];
sx q[3];
rz(12.5771424531858) q[3];
cx q[3],q[2];
rz(1.89697408676147) q[2];
sx q[2];
rz(4.97327950795228) q[2];
sx q[2];
rz(9.66773816048309) q[2];
rz(2.77295589447021) q[3];
sx q[3];
rz(3.74710032542283) q[3];
sx q[3];
rz(8.19254086016818) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.315455615520477) q[0];
sx q[0];
rz(6.06050530274446) q[0];
sx q[0];
rz(9.05745560526058) q[0];
rz(0.796337246894836) q[1];
sx q[1];
rz(2.08347359498078) q[1];
sx q[1];
rz(10.0672815799634) q[1];
cx q[1],q[0];
rz(-1.31430196762085) q[0];
sx q[0];
rz(1.99019423325593) q[0];
sx q[0];
rz(10.4767812251966) q[0];
rz(0.965131402015686) q[2];
sx q[2];
rz(5.64095774491365) q[2];
sx q[2];
rz(10.0118962287824) q[2];
cx q[2],q[1];
rz(0.312562346458435) q[1];
sx q[1];
rz(4.15265837510163) q[1];
sx q[1];
rz(12.0316238164823) q[1];
rz(0.238249972462654) q[3];
sx q[3];
rz(3.63519421418244) q[3];
sx q[3];
rz(10.1988368391912) q[3];
cx q[3],q[2];
rz(-0.639288604259491) q[2];
sx q[2];
rz(4.07720419962937) q[2];
sx q[2];
rz(11.1960013866346) q[2];
rz(-0.227451995015144) q[3];
sx q[3];
rz(4.39122393925721) q[3];
sx q[3];
rz(8.23320708274051) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.745339035987854) q[0];
sx q[0];
rz(2.96629128058488) q[0];
sx q[0];
rz(10.6229496955793) q[0];
rz(-2.06129693984985) q[1];
sx q[1];
rz(3.35586841602857) q[1];
sx q[1];
rz(12.6612436532895) q[1];
cx q[1],q[0];
rz(-1.08523297309875) q[0];
sx q[0];
rz(4.17454627354676) q[0];
sx q[0];
rz(10.6945940017621) q[0];
rz(0.509682238101959) q[2];
sx q[2];
rz(5.55699268181855) q[2];
sx q[2];
rz(10.3496135234754) q[2];
cx q[2],q[1];
rz(0.669029235839844) q[1];
sx q[1];
rz(4.67147103150422) q[1];
sx q[1];
rz(9.96155074834033) q[1];
rz(-0.613118350505829) q[3];
sx q[3];
rz(3.60046753485734) q[3];
sx q[3];
rz(11.0909079074781) q[3];
cx q[3],q[2];
rz(1.06928110122681) q[2];
sx q[2];
rz(2.53719642956788) q[2];
sx q[2];
rz(9.01830855607196) q[2];
rz(1.96299242973328) q[3];
sx q[3];
rz(4.84296897252137) q[3];
sx q[3];
rz(10.5875501394193) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.562263190746307) q[0];
sx q[0];
rz(3.00709061523015) q[0];
sx q[0];
rz(9.08877528309032) q[0];
rz(-0.520403683185577) q[1];
sx q[1];
rz(4.00576737721498) q[1];
sx q[1];
rz(9.32049613296195) q[1];
cx q[1],q[0];
rz(1.26622426509857) q[0];
sx q[0];
rz(4.38384559948976) q[0];
sx q[0];
rz(9.79233331083461) q[0];
rz(-0.0141167975962162) q[2];
sx q[2];
rz(4.5559288581186) q[2];
sx q[2];
rz(10.3845425605695) q[2];
cx q[2],q[1];
rz(0.194600895047188) q[1];
sx q[1];
rz(2.23339882691438) q[1];
sx q[1];
rz(10.5228798150937) q[1];
rz(1.33332288265228) q[3];
sx q[3];
rz(2.67696115572984) q[3];
sx q[3];
rz(10.4016288280408) q[3];
cx q[3],q[2];
rz(2.52295541763306) q[2];
sx q[2];
rz(4.239819200831) q[2];
sx q[2];
rz(7.42773602008029) q[2];
rz(-0.54730224609375) q[3];
sx q[3];
rz(4.36989703972871) q[3];
sx q[3];
rz(8.33714900015994) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.889214098453522) q[0];
sx q[0];
rz(2.95048983593518) q[0];
sx q[0];
rz(9.97915331124469) q[0];
rz(4.05281448364258) q[1];
sx q[1];
rz(1.26340416272218) q[1];
sx q[1];
rz(9.12336344122096) q[1];
cx q[1],q[0];
rz(0.754795551300049) q[0];
sx q[0];
rz(3.78129223187501) q[0];
sx q[0];
rz(9.81366536616489) q[0];
rz(1.18251717090607) q[2];
sx q[2];
rz(4.51975324948365) q[2];
sx q[2];
rz(11.3905837297361) q[2];
cx q[2],q[1];
rz(0.44228783249855) q[1];
sx q[1];
rz(4.5000400861078) q[1];
sx q[1];
rz(9.81389868854686) q[1];
rz(-1.25701117515564) q[3];
sx q[3];
rz(1.21263924439485) q[3];
sx q[3];
rz(10.8977762222211) q[3];
cx q[3],q[2];
rz(0.0214127954095602) q[2];
sx q[2];
rz(3.32259039779241) q[2];
sx q[2];
rz(9.43170079625353) q[2];
rz(0.931124687194824) q[3];
sx q[3];
rz(5.15187898476655) q[3];
sx q[3];
rz(8.31560764311954) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.508650064468384) q[0];
sx q[0];
rz(3.9786361177736) q[0];
sx q[0];
rz(10.9585705757062) q[0];
rz(-0.70266979932785) q[1];
sx q[1];
rz(3.67280581791932) q[1];
sx q[1];
rz(9.72697356938525) q[1];
cx q[1],q[0];
rz(2.58364987373352) q[0];
sx q[0];
rz(1.46902254422242) q[0];
sx q[0];
rz(10.0864681363027) q[0];
rz(-2.78933310508728) q[2];
sx q[2];
rz(4.0279617031389) q[2];
sx q[2];
rz(15.1283144712369) q[2];
cx q[2],q[1];
rz(-1.85408782958984) q[1];
sx q[1];
rz(2.49374279578263) q[1];
sx q[1];
rz(12.2767271757047) q[1];
rz(-0.328656584024429) q[3];
sx q[3];
rz(0.160200031595775) q[3];
sx q[3];
rz(10.1328398346822) q[3];
cx q[3],q[2];
rz(-2.24937844276428) q[2];
sx q[2];
rz(5.23287072976167) q[2];
sx q[2];
rz(10.3508026957433) q[2];
rz(1.21464359760284) q[3];
sx q[3];
rz(3.63480696280534) q[3];
sx q[3];
rz(9.14848220943614) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.723254263401031) q[0];
sx q[0];
rz(5.51619449456269) q[0];
sx q[0];
rz(10.5333008527677) q[0];
rz(0.333798944950104) q[1];
sx q[1];
rz(4.46828666527803) q[1];
sx q[1];
rz(10.5105846881787) q[1];
cx q[1],q[0];
rz(-1.38768243789673) q[0];
sx q[0];
rz(3.36174504657323) q[0];
sx q[0];
rz(11.3673325538556) q[0];
rz(0.407073378562927) q[2];
sx q[2];
rz(2.70866394241387) q[2];
sx q[2];
rz(8.29322514533206) q[2];
cx q[2],q[1];
rz(1.9232417345047) q[1];
sx q[1];
rz(-0.783037988347463) q[1];
sx q[1];
rz(8.65171680449649) q[1];
rz(0.636607527732849) q[3];
sx q[3];
rz(1.45983401139314) q[3];
sx q[3];
rz(9.10889420508548) q[3];
cx q[3],q[2];
rz(2.62560248374939) q[2];
sx q[2];
rz(3.50608536799485) q[2];
sx q[2];
rz(8.26052055358096) q[2];
rz(0.268167555332184) q[3];
sx q[3];
rz(5.04029396374757) q[3];
sx q[3];
rz(10.1825944542806) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.537157297134399) q[0];
sx q[0];
rz(3.66138920386369) q[0];
sx q[0];
rz(11.4173937797467) q[0];
rz(-0.29414290189743) q[1];
sx q[1];
rz(4.70006838639314) q[1];
sx q[1];
rz(10.4672746419828) q[1];
cx q[1],q[0];
rz(0.718716859817505) q[0];
sx q[0];
rz(5.04246798356111) q[0];
sx q[0];
rz(8.65481845139667) q[0];
rz(1.79965686798096) q[2];
sx q[2];
rz(2.3386152108484) q[2];
sx q[2];
rz(5.36809442042514) q[2];
cx q[2],q[1];
rz(-0.52400541305542) q[1];
sx q[1];
rz(4.53890660603578) q[1];
sx q[1];
rz(8.22530159949466) q[1];
rz(0.036829736083746) q[3];
sx q[3];
rz(5.61884442170198) q[3];
sx q[3];
rz(9.43263737614408) q[3];
cx q[3],q[2];
rz(0.597575485706329) q[2];
sx q[2];
rz(5.61210265954072) q[2];
sx q[2];
rz(8.06006572245761) q[2];
rz(-0.652586042881012) q[3];
sx q[3];
rz(3.93288788397843) q[3];
sx q[3];
rz(10.0731317758481) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.5156010389328) q[0];
sx q[0];
rz(2.76432717044885) q[0];
sx q[0];
rz(9.40743299051329) q[0];
rz(0.0167724397033453) q[1];
sx q[1];
rz(2.3544631918245) q[1];
sx q[1];
rz(9.27089249192878) q[1];
cx q[1],q[0];
rz(3.23249244689941) q[0];
sx q[0];
rz(2.81503009994561) q[0];
sx q[0];
rz(11.2332392692487) q[0];
rz(0.88123494386673) q[2];
sx q[2];
rz(3.98526528676087) q[2];
sx q[2];
rz(8.89494094847842) q[2];
cx q[2],q[1];
rz(-0.202752336859703) q[1];
sx q[1];
rz(1.86492636998231) q[1];
sx q[1];
rz(10.6159472227018) q[1];
rz(-0.0758991688489914) q[3];
sx q[3];
rz(1.59469774563844) q[3];
sx q[3];
rz(10.5430193900983) q[3];
cx q[3],q[2];
rz(1.95355701446533) q[2];
sx q[2];
rz(2.32536581357057) q[2];
sx q[2];
rz(8.980271524183) q[2];
rz(-0.19017530977726) q[3];
sx q[3];
rz(4.72515788872773) q[3];
sx q[3];
rz(10.7975193023603) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(1.0012149810791) q[0];
sx q[0];
rz(4.39113322098786) q[0];
sx q[0];
rz(10.4954306840818) q[0];
rz(-0.0459146872162819) q[1];
sx q[1];
rz(4.81179061730439) q[1];
sx q[1];
rz(10.2158483028333) q[1];
cx q[1],q[0];
rz(0.850020289421082) q[0];
sx q[0];
rz(3.48872426350648) q[0];
sx q[0];
rz(9.16471994518443) q[0];
rz(1.17317986488342) q[2];
sx q[2];
rz(4.96102717717225) q[2];
sx q[2];
rz(10.8472896575849) q[2];
cx q[2],q[1];
rz(0.485342711210251) q[1];
sx q[1];
rz(3.68687764008576) q[1];
sx q[1];
rz(8.59028342961475) q[1];
rz(0.505133807659149) q[3];
sx q[3];
rz(4.82871690590913) q[3];
sx q[3];
rz(11.3671284675519) q[3];
cx q[3],q[2];
rz(-1.64482402801514) q[2];
sx q[2];
rz(6.09858551819856) q[2];
sx q[2];
rz(11.0936418533246) q[2];
rz(1.52769267559052) q[3];
sx q[3];
rz(2.08526411850984) q[3];
sx q[3];
rz(10.6649679899137) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.304626733064651) q[0];
sx q[0];
rz(4.00154748757417) q[0];
sx q[0];
rz(9.52084804921552) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(0.216653600335121) q[1];
sx q[1];
rz(6.09936467011506) q[1];
sx q[1];
rz(8.57584199904605) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(2.15474390983582) q[2];
sx q[2];
rz(4.14775231679017) q[2];
sx q[2];
rz(8.83213130234882) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(1.26441061496735) q[3];
sx q[3];
rz(3.77392903168733) q[3];
sx q[3];
rz(8.87922052144214) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
