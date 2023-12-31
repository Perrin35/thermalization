OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.176856696605682) q[0];
sx q[0];
rz(5.4194873889261) q[0];
sx q[0];
rz(9.62588622271224) q[0];
rz(1.74452984333038) q[1];
sx q[1];
rz(4.94645241101319) q[1];
sx q[1];
rz(10.0516047835271) q[1];
cx q[1],q[0];
rz(-0.62030827999115) q[0];
sx q[0];
rz(4.59339454968507) q[0];
sx q[0];
rz(12.7557272672574) q[0];
rz(-1.21885359287262) q[2];
sx q[2];
rz(4.94201138814027) q[2];
sx q[2];
rz(13.2923540830533) q[2];
cx q[2],q[1];
rz(0.991702914237976) q[1];
sx q[1];
rz(3.82580629189546) q[1];
sx q[1];
rz(12.0646514654081) q[1];
rz(0.563431918621063) q[3];
sx q[3];
rz(3.46346655686433) q[3];
sx q[3];
rz(9.67540792226001) q[3];
cx q[3],q[2];
rz(0.255192488431931) q[2];
sx q[2];
rz(4.07237670023973) q[2];
sx q[2];
rz(8.13126287459537) q[2];
rz(-1.10399293899536) q[3];
sx q[3];
rz(3.98934266169602) q[3];
sx q[3];
rz(10.1795713663022) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.247938081622124) q[0];
sx q[0];
rz(5.22413698037202) q[0];
sx q[0];
rz(8.43959692715808) q[0];
rz(2.74241161346436) q[1];
sx q[1];
rz(1.87894109089906) q[1];
sx q[1];
rz(9.31741499005958) q[1];
cx q[1],q[0];
rz(-1.02213144302368) q[0];
sx q[0];
rz(3.8532926758104) q[0];
sx q[0];
rz(8.77109614609882) q[0];
rz(1.15468800067902) q[2];
sx q[2];
rz(2.51187178690965) q[2];
sx q[2];
rz(8.47572711705371) q[2];
cx q[2],q[1];
rz(-0.799223601818085) q[1];
sx q[1];
rz(4.07201269467408) q[1];
sx q[1];
rz(11.4625847101133) q[1];
rz(1.18223416805267) q[3];
sx q[3];
rz(4.86826852162416) q[3];
sx q[3];
rz(10.0985609650533) q[3];
cx q[3],q[2];
rz(2.51695156097412) q[2];
sx q[2];
rz(4.89415696461732) q[2];
sx q[2];
rz(8.76307860612079) q[2];
rz(0.0557569041848183) q[3];
sx q[3];
rz(3.5006019790941) q[3];
sx q[3];
rz(6.52903339862033) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.677028834819794) q[0];
sx q[0];
rz(3.69548800786073) q[0];
sx q[0];
rz(9.46512782051369) q[0];
rz(0.579987466335297) q[1];
sx q[1];
rz(2.2433387358957) q[1];
sx q[1];
rz(10.5071542024533) q[1];
cx q[1],q[0];
rz(1.49899351596832) q[0];
sx q[0];
rz(5.07855454285676) q[0];
sx q[0];
rz(11.5470077752988) q[0];
rz(3.53353404998779) q[2];
sx q[2];
rz(5.90225568612153) q[2];
sx q[2];
rz(8.71188161372348) q[2];
cx q[2],q[1];
rz(0.129900023341179) q[1];
sx q[1];
rz(4.49082055886323) q[1];
sx q[1];
rz(12.8107654809873) q[1];
rz(-1.79369223117828) q[3];
sx q[3];
rz(1.88670233090455) q[3];
sx q[3];
rz(13.5651874303739) q[3];
cx q[3],q[2];
rz(3.00641560554504) q[2];
sx q[2];
rz(5.02949574788148) q[2];
sx q[2];
rz(8.8708853483121) q[2];
rz(1.64843773841858) q[3];
sx q[3];
rz(2.08863833745057) q[3];
sx q[3];
rz(12.013856625549) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.647881388664246) q[0];
sx q[0];
rz(0.584090622263499) q[0];
sx q[0];
rz(9.38939569740697) q[0];
rz(2.14540505409241) q[1];
sx q[1];
rz(4.67459061940248) q[1];
sx q[1];
rz(13.0544659852903) q[1];
cx q[1],q[0];
rz(-0.490220844745636) q[0];
sx q[0];
rz(4.24804869492585) q[0];
sx q[0];
rz(11.0824036359708) q[0];
rz(2.14962244033813) q[2];
sx q[2];
rz(4.20131400425965) q[2];
sx q[2];
rz(7.79349646567508) q[2];
cx q[2],q[1];
rz(1.07798731327057) q[1];
sx q[1];
rz(1.69971528847749) q[1];
sx q[1];
rz(9.26880765556499) q[1];
rz(2.15321254730225) q[3];
sx q[3];
rz(5.12025943596894) q[3];
sx q[3];
rz(7.53851673602268) q[3];
cx q[3],q[2];
rz(2.78672432899475) q[2];
sx q[2];
rz(7.29831138451631) q[2];
sx q[2];
rz(9.93224689959689) q[2];
rz(-1.18217015266418) q[3];
sx q[3];
rz(1.29272619088227) q[3];
sx q[3];
rz(10.6797622203748) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.599471569061279) q[0];
sx q[0];
rz(0.713897140818187) q[0];
sx q[0];
rz(9.05451456307575) q[0];
rz(0.263693571090698) q[1];
sx q[1];
rz(4.71951690514619) q[1];
sx q[1];
rz(9.62129657565757) q[1];
cx q[1],q[0];
rz(-0.580134332180023) q[0];
sx q[0];
rz(3.55973130662973) q[0];
sx q[0];
rz(11.7446293592374) q[0];
rz(-1.60280847549438) q[2];
sx q[2];
rz(4.25150922139222) q[2];
sx q[2];
rz(9.08930111526653) q[2];
cx q[2],q[1];
rz(-2.64247512817383) q[1];
sx q[1];
rz(3.45617729623849) q[1];
sx q[1];
rz(16.6794743299405) q[1];
rz(1.29451775550842) q[3];
sx q[3];
rz(4.87607148488099) q[3];
sx q[3];
rz(9.02921161650821) q[3];
cx q[3],q[2];
rz(-2.12815165519714) q[2];
sx q[2];
rz(5.33148923714692) q[2];
sx q[2];
rz(11.6492721795957) q[2];
rz(0.83459484577179) q[3];
sx q[3];
rz(4.97384265263612) q[3];
sx q[3];
rz(9.09109229444667) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.82814002037048) q[0];
sx q[0];
rz(4.84172180493409) q[0];
sx q[0];
rz(9.22274335323974) q[0];
rz(2.06878852844238) q[1];
sx q[1];
rz(4.64830425580079) q[1];
sx q[1];
rz(4.53540942668124) q[1];
cx q[1],q[0];
rz(2.41123390197754) q[0];
sx q[0];
rz(3.48200175364549) q[0];
sx q[0];
rz(11.2254429817121) q[0];
rz(-0.27758526802063) q[2];
sx q[2];
rz(3.52958387334878) q[2];
sx q[2];
rz(9.27854745685264) q[2];
cx q[2],q[1];
rz(1.26536333560944) q[1];
sx q[1];
rz(2.66370421846444) q[1];
sx q[1];
rz(8.81944284438297) q[1];
rz(0.835157811641693) q[3];
sx q[3];
rz(5.0613501389795) q[3];
sx q[3];
rz(8.84349594115421) q[3];
cx q[3],q[2];
rz(0.949719250202179) q[2];
sx q[2];
rz(4.62159517605836) q[2];
sx q[2];
rz(8.58246246575519) q[2];
rz(-1.08742821216583) q[3];
sx q[3];
rz(2.4098586161905) q[3];
sx q[3];
rz(11.0833542108457) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.61206459999084) q[0];
sx q[0];
rz(5.01351419289643) q[0];
sx q[0];
rz(11.5450088739316) q[0];
rz(1.13133239746094) q[1];
sx q[1];
rz(5.33349839051301) q[1];
sx q[1];
rz(12.7766217946927) q[1];
cx q[1],q[0];
rz(0.700498580932617) q[0];
sx q[0];
rz(2.99196176429326) q[0];
sx q[0];
rz(9.88470423816844) q[0];
rz(-0.0657966732978821) q[2];
sx q[2];
rz(4.87571862538392) q[2];
sx q[2];
rz(10.2171407103459) q[2];
cx q[2],q[1];
rz(-0.103403128683567) q[1];
sx q[1];
rz(1.25451711018617) q[1];
sx q[1];
rz(6.08172771929904) q[1];
rz(0.639560580253601) q[3];
sx q[3];
rz(2.31244716246659) q[3];
sx q[3];
rz(11.6862325429837) q[3];
cx q[3],q[2];
rz(-2.31470990180969) q[2];
sx q[2];
rz(4.89333167870576) q[2];
sx q[2];
rz(10.805917954437) q[2];
rz(0.7621049284935) q[3];
sx q[3];
rz(2.66571811039979) q[3];
sx q[3];
rz(9.83823878168269) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.49522265791893) q[0];
sx q[0];
rz(0.782839926081248) q[0];
sx q[0];
rz(8.83753625153705) q[0];
rz(-3.5878632068634) q[1];
sx q[1];
rz(4.4210286458307) q[1];
sx q[1];
rz(8.4480577468793) q[1];
cx q[1],q[0];
rz(1.85510885715485) q[0];
sx q[0];
rz(3.62033906777436) q[0];
sx q[0];
rz(11.3229807376783) q[0];
rz(0.459372103214264) q[2];
sx q[2];
rz(2.96152353485162) q[2];
sx q[2];
rz(10.5510475397031) q[2];
cx q[2],q[1];
rz(-2.63736200332642) q[1];
sx q[1];
rz(5.14978018601472) q[1];
sx q[1];
rz(10.6113358497541) q[1];
rz(1.01838898658752) q[3];
sx q[3];
rz(4.07980677683885) q[3];
sx q[3];
rz(7.75007555483981) q[3];
cx q[3],q[2];
rz(-1.74051594734192) q[2];
sx q[2];
rz(3.84986433585221) q[2];
sx q[2];
rz(9.97722074984714) q[2];
rz(-0.465941220521927) q[3];
sx q[3];
rz(4.8955067714029) q[3];
sx q[3];
rz(11.5388595819394) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.47336459159851) q[0];
sx q[0];
rz(2.34782740672166) q[0];
sx q[0];
rz(7.20383474826022) q[0];
rz(3.09055113792419) q[1];
sx q[1];
rz(4.69811716874177) q[1];
sx q[1];
rz(8.52568260430499) q[1];
cx q[1],q[0];
rz(0.162227869033813) q[0];
sx q[0];
rz(5.5290962775522) q[0];
sx q[0];
rz(10.0698747992437) q[0];
rz(-2.58286833763123) q[2];
sx q[2];
rz(4.99710288842256) q[2];
sx q[2];
rz(8.77189437150165) q[2];
cx q[2],q[1];
rz(-1.85832619667053) q[1];
sx q[1];
rz(2.79106334050233) q[1];
sx q[1];
rz(9.08387140034839) q[1];
rz(1.57100975513458) q[3];
sx q[3];
rz(5.08256140549714) q[3];
sx q[3];
rz(9.82571820019885) q[3];
cx q[3],q[2];
rz(0.945453226566315) q[2];
sx q[2];
rz(3.07101994951303) q[2];
sx q[2];
rz(9.32018922119542) q[2];
rz(0.83834046125412) q[3];
sx q[3];
rz(3.89090827305848) q[3];
sx q[3];
rz(11.683614230148) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.339124411344528) q[0];
sx q[0];
rz(3.95925733645494) q[0];
sx q[0];
rz(10.5426060914914) q[0];
rz(0.712401866912842) q[1];
sx q[1];
rz(5.65196314652497) q[1];
sx q[1];
rz(9.82174268960162) q[1];
cx q[1],q[0];
rz(-0.944969415664673) q[0];
sx q[0];
rz(5.19016757805879) q[0];
sx q[0];
rz(9.96899483203098) q[0];
rz(2.77515411376953) q[2];
sx q[2];
rz(1.28411844571168) q[2];
sx q[2];
rz(8.5888978600423) q[2];
cx q[2],q[1];
rz(1.1399998664856) q[1];
sx q[1];
rz(4.62401023705537) q[1];
sx q[1];
rz(9.08786672949001) q[1];
rz(2.47305369377136) q[3];
sx q[3];
rz(6.88881483872468) q[3];
sx q[3];
rz(12.2521419286649) q[3];
cx q[3],q[2];
rz(4.08999967575073) q[2];
sx q[2];
rz(2.24141749938066) q[2];
sx q[2];
rz(6.18653104304477) q[2];
rz(1.72766733169556) q[3];
sx q[3];
rz(5.14513388474519) q[3];
sx q[3];
rz(11.4394280671994) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.142508983612061) q[0];
sx q[0];
rz(4.96524682839448) q[0];
sx q[0];
rz(9.89725646971866) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(3.23650574684143) q[1];
sx q[1];
rz(3.94865003426606) q[1];
sx q[1];
rz(6.57678053378269) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(2.96988248825073) q[2];
sx q[2];
rz(3.43815526564653) q[2];
sx q[2];
rz(6.57487199305698) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.602091073989868) q[3];
sx q[3];
rz(2.76293039520318) q[3];
sx q[3];
rz(8.66674069165393) q[3];
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
