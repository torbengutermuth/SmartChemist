"""smart_chemist URL definitions"""
from django.urls import path
from rest_framework.routers import DefaultRouter
from . import views

urlpatterns = [
    path('api/names', views.SmartChemistView.as_view()),
    path('api/bugs', views.SmartChemistBugReportView.as_view()),
]

router = DefaultRouter()
router.register('jobs', views.PatternMatchingJobViewSet)
router.register('jobs_input', views.PatternMatchingInputModelViewSet)
router.register('jobs_output', views.PatternMatchingOutputModelViewSet)
urlpatterns.extend(router.urls)
